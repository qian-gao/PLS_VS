# PLS_VS

PLSvar_sel aims to provide a variable selection tool for data derived from metabolomic studies with a time-series design. It models the data based on multilinear partial least squares regression (NPLS) with different combinations of bilinear/trilinear X and group/time-response dummy Y. A bootstrap procedure is then performed to calculate the average VIP for each variable forming the basis for variable selection.

This code was developed based on the NPLS package (Bro, R. Multiway calibration. Multilinear PLS. J. Chemom. 1996, 10, 47–61) and the VIP_nway package (Favilla, S.; Durante, C.; Vigni, M. L.; Cocchi, M. Assessing feature relevance in NPLS models by VIP. Chemom. Intell. Lab. Syst. 2013, 129, 76–86).
