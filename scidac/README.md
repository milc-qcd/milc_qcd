### Compile

For gb_baryon applications, only QMP and QIO are needed. Fisrt, you need to change the configuration files at ```qinstall/qio```
and ```qinstall/qmp```, especially specifying the C compiler with the ```CC``` option.

Once this is done, you can access the folder ```qio``` and do

```sh
make install_qio
```
 which will install both QMP and QIO to the ```install``` folder (always with the suffix ```-openmpi``` but it doesn't mean
 it is compiled with openmpi).
