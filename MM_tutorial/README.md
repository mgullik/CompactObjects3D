Here you'll find power, cross, lags spectra and coherence function for a number of sources and missions with the corresponding .xcm files to load the files and the models in Xspec. Each subdirectlry has a README.md file with the description of the files available.

**Start with the directory [RXTE_GRS_1915+105](https://github.com/EdNathan/CompactObjects3D/tree/main/MM_tutorial/RXTE_GRS_1915%2B105), since the .xcm files have many comments explaining what it is done.**

In general the directories contain:

Power spectrum band 1 (PS1)
Power spectrum band 2 (PS2)
Real part of the cross spectrum of band 1 vs. band 2 (RE)
Imaginary part of the cross spectrum of band 1 vs. band 2 (IM)
Phase lags of band 1 vs. band 2 (PL)
Coherence function of band 1 vs. band 2 (CF)
The XSPEC .xcm file(s) to read the data and define the model 

Normally you would go to the directpry with all the files, run XSPEC and then run the .xcm file:

```
\> xspec
XSPEC12> @file.xcm
```

