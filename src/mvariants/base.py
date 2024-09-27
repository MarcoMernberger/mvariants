#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""base.py: Contains some base modules needed in caller.py and preprocessor.py."""

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any
from mvariants import __version__
from mbf.externals.util import download_zip_and_turn_into_tar_gzip
import warnings
import subprocess
import re


class OptionHandler:
    """
    Parameter handling class that performs basic checks for given paremeters.

    Convenience class that ensures some basic option handling methods for
    external algorithms, including a simple parameter checking. Must be
    initalized with a function that returns allowed and default
    parameters for an external method, as this is specific for the method
    to be used.

    Parameters
    ----------
    name : str
        Unique id of the instance.
    option_getter : Callable[[], Tuple[str, Dict, Dict]]
        Callable that returns the help string, a dictionary of allowed options
        and descriptions and a dictionary with default values.
    warn : st, optional
        Determines how the option handler treats incorrect parameters, by
        default "always". This value is used as warnings filter option.
        See python package warning for details.
    explicit : bool, optional
        Determines if default parameters inferred from the help string should
        be set explicitly, by default False.
    """

    def __init__(
        self,
        name: str,
        option_getter: Callable[[], Tuple[str, Dict, Dict]],
        warn: str = "always",
        explicit: bool = False,
    ) -> None:
        """OptionHandler constructor, see class documentation for details."""
        self.name = name
        self.warn = warn
        self.option_getter = option_getter
        self.explicit = explicit
        self.__initialize_dicts()

    def __initialize_dicts(self) -> None:
        """Sets internal attributes by calling the option getter function."""
        (
            self._help_str,
            self._accepted_arguments,
            self._default_parameter,
        ) = self.option_getter()

    @property
    def default_parameter(self) -> Dict[str, Any]:
        """
        Getter for the default parameters of the external methods.

        Returns a dictionary containing default parameters and values for an
        external method, initially returned by the option_getter.

        Returns
        -------
        Dict[str, Any]
            Dictionary of default parameters and values.
        """
        return self._default_parameter

    @default_parameter.setter
    def default_parameter(self, options: Dict[str, Any]) -> None:
        """Setter for the default parameters."""
        self._default_parameter = options

    @property
    def accepted_arguments(self) -> Dict[str, str]:
        """
        Returns allowed parameters for an external method.

        Returns a dictionary of allowed options for an external method and
        their corresponding descriptions.

        Returns
        -------
        Dict[str, str]
            Dictionary of allowed options and corresponding decriptions.
        """
        return self._accepted_arguments

    @accepted_arguments.setter
    def accepted_arguments(self, options: Dict[str, str]):
        """Setter for accepted arguments."""
        self._accepted_arguments = options

    def accepted_arguments_str(self) -> str:
        """
        Returns allowed options as a string.

        Returns a string representation of the accepted_arguments attribute.

        Returns
        -------
        str
            String of accepted arguments.
        """
        return "\n".join([f"{item[0]}:{item[1]}" for item in self.accepted_arguments.items()])

    def print_help(self) -> None:
        """Prints the help string of the external method."""
        print(self._help_str)

    def check_options(self, options: Dict[str, Any]) -> Dict[str, Any]:
        """
        Performs a basic check of parameters supplied to an external method.

        Checks a dictionary of supplied parameters for an external method and
        compares it to accepted arguments. If set to explicit, this will
        also add unspecified parameters to the options with default settings.
        This does not check the supplied values itself, only the parameter
        names. By default, unknown parameters will be removed and a Warning is
        raised.

        Parameters
        ----------
        options : Dict[str, Any]
            List of parameter for the external method.

        Returns
        -------
        Dict[str, Any]
            Updated list of parameters.

        Raises
        ------
        UserWarning
            User warning for unknown parameters in the options.
        """
        warnings.simplefilter(self.warn, UserWarning)
        if self.accepted_arguments is None:
            warnings.warn(
                f"{self.name}: option_getter returned no accepted_arguments, no option checking is performed."
            )
            return options
        else:
            if self.explicit:
                checked = self.default_parameter
            else:
                checked = {}
            for arg in options:
                if arg not in self._accepted_arguments:
                    warnings.warn(
                        f"{arg} is not a valid parameter, will be ignored. Use any of: \n{list(self.accepted_arguments.keys())}."
                    )
                else:
                    checked[arg] = options[arg]
            return checked

    @staticmethod
    def options_as_list(options: Dict[str, Any]) -> List[Any]:
        """
        Turns a dictionary of parameters as a list.

        Parameters
        ----------
        options : Dict[str, Any]
            Dictionary with parameters.

        Returns
        -------
        List[Any]
            Flattened list of parameters and values.
        """
        ret = []
        for pair in options.items():
            if pair[1] == "":
                ret.append(pair[0])
            elif isinstance(pair[1], list):
                for x in pair[1]:
                    ret.extend([pair[0], x])
            else:
                ret.extend([pair[0], pair[1]])
        return ret

    @staticmethod
    def options_as_list_str(options: Dict[str, Any]) -> List[str]:
        """
        Returns options as a list of strings.

        Flattens a supplied dictionary with parameters and values to a list
        and turns it into string representations to be supplied to subprocess.

        Parameters
        ----------
        options : Dict[str, Any]
            Dictionary of parameters and values.

        Returns
        -------
        List[str]
            Flattened list of parameters and values as strings.
        """
        return [str(x) for x in OptionHandler.options_as_list(options)]


class GATK:
    """
    Wrapper class for the GATK toolbox.

    Wrapper for the GATK toolbox that subclasses mbf.externals.ExternalAlgorithm
    to take care of version handling.
    Allows to invoke any GATK method and return help strings and option
    handling.

    Parameters
    ----------
    tool : Optional[str], optional
        The name of the GATK tool to be invoked, by default None. Set to None
        if either none or multiple tools will be used. If None, options are
        inferred from the basic command call.
    options : Optional[Dict[str, str]], optional
        Dictionary of parameter-values to be supplied to the GATK call.
    version : str, optional
        Version of the tool to be used, by default "_last_used".
    store : ExterrnalAlgorithmStore, optional
        Store that handles downloads and provides thh external tool, by default
        None.

    Raises
    ------
    ValueError
        Raises if the supplied GATK command is not valid.
    """

    def __init__(
        self,
        tool: Optional[str] = None,
        options: Optional[Dict[str, str]] = {},
        **kwargs,
    ):
        """GATK constructor, see class documentation for details."""
        defname = "base"
        self.command = "gatk"
        tools = self.parse_tools()
        if tool is not None:
            if tool not in tools:
                raise ValueError(f"Unknown tool {tool}, use 'None' or any of {tools}.")
            defname = tool
        self.tool = tool
        self.instance_name = kwargs.get("instance_name", "_".join(["GATK", defname]))
        explicit = kwargs.get("explicit_options", True)
        # self.optionhandler = OptionHandler(
        #    f"{self.instance_name}_optionhandler",
        #    self.get_option_parser(),
        #    warn="always",
        #    explicit=explicit,
        # )
        # self.options = self.optionhandler.check_options(options)
        # if len(self.options) == 0:
        #    self.options = self.optionhandler.default_parameter
        self.options = options

    def fetch_version(self, version: str, target_filename: Path) -> None:
        """
        Takes care of the tool download.

        Overrides the ExternalAlgorithm methood. Downloads the external method
        to the prebuild location specified by the corresponding
        ExternalAlgorithmStore and packs it into a tar.gz file.

        Parameters
        ----------
        version : str
            The tool version to be used.
        target_filename : Path
            The path to the local tar.gz file.
        """
        url = (
            f"https://github.com/broadinstitute/gatk/releases/download/{version}/gatk-{version}.zip"
        )
        download_zip_and_turn_into_tar_gzip(
            url, target_filename, chmod_x_files=[f"gatk-{version}/gatk"]
        )

    @property
    def name(self) -> str:
        """
        Returns the name of the external method for version handling.

        Overrides the ExternalAlgorithm method.

        Returns
        -------
        str
            Name of the external method.
        """
        return "GATK"

    def print_help(self) -> None:
        """
        Prints a help string for the external method.

        Calls the OptionHandler function to print a help string of the
        corresponding OPtionHandler object.
        """
        self.optionhandler.print_help()

    def print_tools(self) -> None:
        """
        Prints a list of accepted GATK tools.

        Calls the GATK tool to get a list of accepted tool commands and prints
        them.

        [extended_summary]
        """
        self.run_command(["--list"])

    @property
    def multi_core(self) -> bool:
        """
        Returns True if the GATK call can use multiple cores.

        Overrides the ExternalAlgorithm method.
        """
        return "-nct" in self.options or "-nt" in self.options

    def build_cmd(self, output_directory: Optional[Path], ncores: int, arguments: List[str]):
        """
        Returns a command as a list of strings to be passed to subprocess.

        Constructs an executable command to be passed to subprocess from the
        output directory, the number of cores to be used and a list of command
        options. Overrides the ExternalAlgorithm method and ignores the output
        directory and ncores parameters. Both are handled by the option list.

        Parameters
        ----------
        output_directory : Optional[Path]
            Path of output directory.
        ncores : int
            Number of cores to be used.
        arguments : List[str]
            List of string arguments for the command call.

        Returns
        -------
        List[str]
            Command to subprocess as a list of strings.
        """
        if self.tool is None:
            return [self.command] + arguments
        return [self.command, self.tool] + arguments

    def __repr__(self) -> str:
        return f"GATK({self.tool}, {self.options})"

    def __str__(self) -> str:
        return f"GATK(tool = {self.tool}, options = {self.options})"

    def get_cores_needed(self) -> int:
        """
        Returns the number of cores the GATK call will need.

        Overrides the ExternalAlgorithm abstract method.

        Returns
        -------
        int
            [description]
        """
        if "-nt" in self.options:
            return self.options["-nct"]
        elif "-nct" in self.options:
            return self.options["-nct"]
        return 1

    def parse_tools(self) -> List[str]:
        """
        Returns a list of accepted GATK tools.

        Calls the GATK toolkit with the --list option and parses the result
        string for allowed tools.

        Returns
        -------
        List[str]
            List of tool options for GATK.
        """
        command = [self.command, "--list"]
        try:
            help_str = subprocess.check_output(
                command,
                stderr=subprocess.STDOUT,
            ).decode()
        except subprocess.CalledProcessError as exc:
            help_str = exc.output.decode()
        tools = []
        for line in help_str.split("\n"):
            if len(line) > 0:
                if line.startswith("\x1b[32m"):
                    rescape = re.compile(r"\x1b[^m]*m")
                    line = rescape.sub("", line)
                    tools.append(line.strip().split()[0])
        return tools

    def get_option_parser(self) -> Callable:
        """
        Returns a callable to be passed to the OptionHandler.

        Returns a callable that returns the help string, a dictionary of
        allowed options and descriptions and a dictionary with default values.

        Returns
        -------
        Callable
            option_getter for the corresponding OptionHandler.
        """

        def __init_options() -> Tuple[str, Dict[str, str], Dict[str, Any]]:
            command = self.build_cmd(None, 1, ["--help"])
            print(command)
            try:
                help_str = subprocess.check_output(
                    command,
                    stderr=subprocess.STDOUT,
                ).decode()
            except subprocess.CalledProcessError as exc:
                help_str = exc.output.decode()
            allowed_options = {}
            default_values = {}
            valid = ""

            def __commit(line, default_values, allowed_options, valid):
                if line is None:
                    pass
                elif line.startswith("--"):
                    splits = re.split(r"\s{1,}", line.strip())
                    print("splits", splits)
                    splits = [splits[0], " ".join(splits[1:] + [valid])]
                    required = "Required." in splits[1]
                    opt = splits[0]
                    sep = opt.find(",")
                    if sep > 0:
                        opt = opt[sep + 1 :]
                        print(opt)
                    option_name, option_type = tuple(opt.split(":"))
                    allowed_options[option_name] = splits[1] + f"Type={option_type}."
                    if not required and not valid.startswith("Valid"):
                        def_index = splits[1].find("Default value: ")
                        if def_index > 0:
                            val = splits[1][def_index + 15 :].split(". ")[0]
                            if val != "null":
                                if val == "":
                                    pass
                                else:
                                    try:
                                        val = eval(val)
                                    except NameError:
                                        pass
                                default_values[option_name] = val
                                if "This argument may be specified 0 or more times" in splits[1]:
                                    default_values[option_name] = list(default_values[option_name])
                else:
                    pass
                return valid

            newline = ""
            for line in help_str.split("\n"):
                print(line)
                if line == "":
                    __commit(newline, default_values, allowed_options, valid)
                    newline = line
                    continue
                elif line.startswith("Valid ") or line.startswith("Advanced Arguments"):
                    valid = line[:-1] + ". "
                else:
                    newline += line
            valid = __commit(newline, default_values, allowed_options, valid)
            if "--version" in default_values:
                del default_values["--version"]
            return help_str, allowed_options, default_values

        if hasattr(self, "_option_parser_modifier"):
            __init_options = self._option_parser_modifier()(__init_options)
        return __init_options

    def run_command(
        self,
        command_arguments: List[str],
        ncores: Optional[int] = 1,
        output_directory: Optional[Path] = None,
    ) -> None:
        """
        Calls GATK toolbox via subprocess.

        Simple command call to GATK via subprocess.

        Parameters
        ----------
        output_directory : Optional[Path]
            Path of output directory.
        ncores : int
            Number of cores to be used.
        arguments : List[str]
            List of string arguments for the command call.
        """
        cmd = self.build_cmd(output_directory, ncores, command_arguments)
        subprocess.check_call(cmd)
