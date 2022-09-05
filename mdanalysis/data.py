import pandas as pd


class Data(pd.DataFrame):

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)