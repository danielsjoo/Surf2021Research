import ctypes

arrayDeltaG = ctypes.CDLL('./c_routines/deltaG.so')

arrayDeltaG.deltaG.restype = None
arrayDeltaG.deltaG.argtypes = (ctypes.c_int, 
                               ctypes.POINTER(ctypes.c_double))

def helper2(sequence):
    #key: 1:A 2:G 3:C 4:T
    x = []
    for char in sequence:
        if char == "A":
            x.append(1)
        elif char == "G":
            x.append(2)
        elif char == "C":
            x.append(3)
        else:
            x.append(4)
    array_size = len(sequence)
    x = (ctypes.c_double * array_size)(*x)
    arrayDeltaG.deltaG(array_size, x)
    return x

def return_sig_positions(sequence):
    x = []
    for i in helper2(sequence)[:]:
        x.append(i)
    return x