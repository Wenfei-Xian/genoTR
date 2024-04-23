# genoTR
A tool for genotyping short tandem repeat.

![alt text](https://github.com/Wenfei-Xian/genoTR/blob/main/Picture1.png)

before install, you should install htslib, spoa and ssw

install
```
g++ genoTR.cpp -o genoTR -Ihtslib/include -Lhtslib/lib -lhts -Ispoa/include -Lspoa/lib -lspoa -lz -Wall ssw_cpp.cpp ssw.c
```
