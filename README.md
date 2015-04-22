Calculates global image statistics by combining sums into groups so as to minimize floating point precision errors

At 3.2 arcsec, the global image is 405000 x 202500 ~ 10^11 pixels 

64-bit double precision floating point number has precision of  15-17 significant decimal digits
32-bit single precision floating point number has precision of 6-9 significant decimal digits:

So, if we are summing over the entire image one pixel at a time, 32-bit float is definitely not enough

double precision may suffice by doing simple sum, but we do a 2-step aggregation here.

for counter data type, we can use long long (64 bit integer) which has range of ~ -10^19 to 10^19
