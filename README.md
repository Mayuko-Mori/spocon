To calculate spot contrast for given spot temperature, photosphere temperature and filters, based on BT-Settl models.
This code is applicable to spot/photosphere temperatures from 1700K to 5000K.

# Usage
```python
import spocon
spocon.contrast(3215, 3350, 'g') #Tspot, Tphot, filter
```
Tspot and Tphot can also be entered in array.

```python
spocon.contrast([3050, 3100, 3150], [3350, 3390], 'g')
```


## Filters
```python
bands = ['g', 'r', 'i', 'z',  #Sinistro
         'zg', 'zr',  #ZTF
         'k', 't', #Kepler, TESS
         'Mg', 'Mr', 'Mi', 'Mz'] #MuSCAT
```
