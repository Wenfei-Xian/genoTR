# genoTR
A tool for genotyping short tandem repeat.

![alt text](https://github.com/Wenfei-Xian/genoTR/blob/main/Picture1.png)

# install
before install, you should install htslib, spoa and ssw

```
g++ genoTR.cpp -o genoTR -Ihtslib/include -Lhtslib/lib -lhts -Ispoa/include -Lspoa/lib -lspoa -lz -Wall ssw_cpp.cpp ssw.c
```

# limitation
This tool is currently only suitable for highly inbred species, as it lacks a phasing stage.
