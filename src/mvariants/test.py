class Test1:

    def __init__(self):
        self.name = "Test1"
        self.value = 2

    def get_wrapper(self):
        def wrapper(func):
            def inner(*args, **kwargs):
                new_args = [args[0] + self.value, args[1] + self.value]
                return func(*new_args, **kwargs)

            return inner

        return wrapper


class Test2:

    def __init__(self, mod_class):
        self.name = "Test2"
        self.mod_class = mod_class

    def return_someting(self):
        def do_something(x, y):
            return x + y

        wrapper = self.mod_class.get_wrapper()
        return wrapper(do_something)


if __name__ == "__main__":

    test1 = Test1()
    test2 = Test2(test1)

    func = test2.return_someting()
    print(func(3, 4))
