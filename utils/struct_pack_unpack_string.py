# Tracer code for packing/unpacking strings to a struct

import struct

def read(format,file):
    """Simplifies reading struct formatted strings from file"""
    return struct.unpack(format,f.read(struct.calcsize(format)))[0]

def write(format,var,file):
    """Simplifies writing to a file using struct compression"""
    file.write(struct.pack(format,var))

with open('out','wb') as f:
    string = "hello"
    write('!b',len(string),f)
    write('!'+str(len(string))+'s',string,f)

with open('out','rb') as f:
    l = read('!b',f)
    print read('!'+str(l)+'s',f)
