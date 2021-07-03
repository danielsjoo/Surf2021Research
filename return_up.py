import ctypes

libsortarray = ctypes.CDLL('./c_routines/libsortarray.so')

libsortarray.sortArray.restype = None
libsortarray.sortArray.argtypes = (ctypes.c_int, 
                                   ctypes.POINTER(ctypes.c_double))

def helper(sequence):
    """Return a sorted copy of the input array or sequence."""
    #key: 1:A 2:G 3:C 4:T
    x = []
    for char in sequence:
        if char == "A":
            x.append(1)
        elif char == "G":
            x.append(2)
        elif char == "C":
            x.append(3)
        elif char == "U":
            x.append(4)
    array_size = len(sequence)
    x = (ctypes.c_double * array_size)(*x)
    libsortarray.sortArray(array_size, x)
    return x

def return_up(sequence):
    x = []
    for i in helper(sequence)[:]:
        x.append(i)
    return x

if __name__ == '__main__':
    seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"
    print ('Unsorted Array:\n', seq)
    print ('Sorted Array:\n', return_up(seq)[:])