# spocon.contrast
To calculate spot contrast for given spot temperature, photosphere temperature and filters, based on BT-Settl models.
This code is applicable to spot/photosphere temperatures from 1700K to 5000K.

## Usage
```python
import spocon
spocon.contrast(3215, 3350, 'g') #Tspot, Tphot, filter
```
Tspot and Tphot can also be entered in array.

```python
spocon.contrast([3050, 3100, 3150], [3350, 3390], 'g')
```

## Available Filters
```python
bands = ['g', 'r', 'i', 'z',  #Sinistro
         'zg', 'zr',  #ZTF
         'k', 't', #Kepler, TESS
         'Mg', 'Mr', 'Mi', 'Mz'] #MuSCAT
```

# spocon.estimate_Tspot_emp
To estimate the spot temperature when you only know the effective temperature of the star.

- model 'Herbst1' : Eq.(4) in Herbst+21, derived from the sample in Berdyugina05
- model 'Herbst2' : Eq.(5) (typo in the paper) in Herbst+21, derived from the sample in Berdyugina05, Biazzo+06, Valio17

## Usage
```python
T_spot_low, T_spot_mid, T_spot_high = spocon.estimate_Tspot_emp(3200, model='Herbst1') #Tphot, model
```

# Notes
More example in "notebooks" directory.

[2024 Sep 13]  To be updated...
- derive modulation amplitude from Tphot, Tspot, fspot, filter
- estimate fspot from Tphot, filter, amplitude
- estimate fspot from Tphot, TESS light curve
- spocon.filter to show filters
- add more explanation when you run spocon.contrast?
- add more filters
- add more stellar models
- update empirical relation to estimate Tspot (Mori et al. in prep)
