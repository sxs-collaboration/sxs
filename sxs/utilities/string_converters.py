
# Assemble some helper functions

def floater(x):
    import numpy as np
    try:
        f = float(x)
    except:
        f = np.nan
    return f

def floaterbound(x):
    import numpy as np
    try:
        f = float(x)
    except:
        try:
            f = float(x.replace("<", ""))
        except:
            f = np.nan
    return f

def norm(x):
    import numpy as np
    try:
        n = np.linalg.norm(x)
    except:
        n = np.nan
    return n

def three_vec(x):
    import numpy as np
    try:
        a = np.array(x, dtype=float)
        if a.shape != (3,):
            raise ValueError("Don't understand input as a three-vector")
    except:
        a = np.array([np.nan, np.nan, np.nan])
    return a

def datetime_from_string(x):
    import pandas as pd
    try:
        dt = pd.to_datetime(x).tz_convert("UTC")
    except:
        dt = pd.to_datetime("1970-1-1").tz_localize("UTC")
    return dt
