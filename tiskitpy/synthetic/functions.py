import numpy as np


def to_DBs(inp):
    """
    Converts values to dBs ref 1

    Args:
        inp (float, list, or np.ndarray): values to convert
    """
    if isinstance(inp, (list, tuple)):
        inp = np.array(inp)
    return 20 * np.log10(inp)


def from_DBs(inp):
    """
    Converts values from dBs ref 1

    Args:
        inp (float, list or np.ndarray): values to convert
    """
    if isinstance(inp, (list, tuple)):
        inp = np.array(inp)
    return np.power(10., inp / 20)


if __name__ == "__main__":
    # Show an example
    print(f'{to_DBs(20000)=}')
    print(f'{from_DBs(25)=}')
