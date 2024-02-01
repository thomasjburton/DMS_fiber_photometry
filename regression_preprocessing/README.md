## Files

*Exp2_Reversal_analysis.py*
<br />
Adapted from [code](https://github.com/crnolan/behapy/blob/main/examples/analyse.py) originally written by Chris Nolan: 
This script was used to analyse pre-processed fiber photometry data organised into the BIDS structure. It should sit in a ""*scripts*" subfolder in the top level of the BIDS structure. e.g.
```
    .\Exp2_Reversal_BIDS\scripts
```

Specifically, this script calculated dFF, handled event mapping and executed the multiple linear regressions on the data from Experiment 2 described [here](https://doi.org/10.1101/2022.01.31.478585) and available [here](https://doi.org/10.6084/m9.figshare.19083647.v2)

Pre-processing of raw fiber photometry data and organisation into the BIDS structure was achived via the [behapy](https://github.com/crnolan/behapy/tree/main) toolkit.



    
