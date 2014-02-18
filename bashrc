#!/bin/bash

export LEMOSEXT=$FOAM_SRC/LEMOS-2.3.x
export LIB_LEMOS_SRC=$FOAM_SRC/LEMOS-2.3.x/libLEMOS-2.3.x
export PATH=$PATH:$LEMOSEXT/bin


alias lemos='cd $LEMOSEXT'
alias lemossol='cd $LEMOSEXT/solver'
alias lemosutil='cd $LEMOSEXT/utils'
