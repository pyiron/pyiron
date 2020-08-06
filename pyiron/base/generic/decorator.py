import inspect 
from pyiron.base.settings.generic import Settings


s = Settings()


def get_default_args(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }


def get_arguments(funct, **kwargs):
    default_arg_dict = get_default_args(func=funct)
    for k in default_arg_dict.keys():
        if k in kwargs.keys():
            try:  # Some objects we can not compare like the atoms structure 
                if kwargs[k] == default_arg_dict[k] and k in s._user_settings.keys():
                    kwargs[k] = s._user_settings[k]
            except AssertionError:  
                pass
        elif k not in kwargs.keys() and k in s._user_settings.keys():
            kwargs[k] = s._user_settings[k]
        elif k not in kwargs.keys() and k not in s._user_settings.keys():
            kwargs[k] = default_arg_dict[k]
    return kwargs


class pyironfunction(object):
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        kwargs = get_arguments(funct=self.f, **kwargs)
        return self.f(*args, **kwargs)
