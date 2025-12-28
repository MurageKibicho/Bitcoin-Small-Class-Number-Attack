#!/usr/bin/env bash

# See http://www.flintlib.org/downloads.html
export FLINT_VERSION=3.4.0

# See https://www.mpfr.org/mpfr-current/#download
export MPFR_VERSION=4.2.2

# See https://www.mpfr.org/mpfr-current/#download
export GMP_VERSION=6.3.0


export BUILD=`pwd`/build
export PREFIX=$BUILD/local

# Create the build directory
mkdir -p $BUILD

set -ev
cd $BUILD
### Download and build GMP
curl https://gmplib.org/download/gmp/gmp-$GMP_VERSION.tar.xz -o gmp-$GMP_VERSION.tar.xz
tar xvf gmp-$GMP_VERSION.tar.xz
cd gmp-$GMP_VERSION
./configure --prefix=$PREFIX
make -j$(nproc)
make install

### Download and build MPFR
cd $BUILD
curl https://www.mpfr.org/mpfr-current/mpfr-$MPFR_VERSION.tar.xz -o mpfr-$MPFR_VERSION.tar.xz
tar xvf mpfr-$MPFR_VERSION.tar.xz
cd mpfr-$MPFR_VERSION
./configure --prefix=$PREFIX --with-gmp=$PREFIX
make -j$(nproc)
make install

### Download and build Flint
cd $BUILD
curl https://flintlib.org/download/flint-$FLINT_VERSION.tar.gz -o flint-$FLINT_VERSION.tar.gz
tar xvf flint-$FLINT_VERSION.tar.gz
cd flint-$FLINT_VERSION
./configure --prefix=$PREFIX --with-gmp=$PREFIX --with-mpfr=$PREFIX
make -j$(nproc)
make install

### Update shell config to set PATH permanently
SHELL_RC=""
if [[ "$SHELL" == */bash ]]; then
    SHELL_RC="$HOME/.bashrc"
elif [[ "$SHELL" == */zsh ]]; then
    SHELL_RC="$HOME/.zshrc"
else
    echo "Unknown shell. Please manually add $PREFIX/bin to your PATH."
fi

if [[ -n "$SHELL_RC" ]]; then
    # Only add the lines if they aren't already present
    grep -qxF "export PATH=\"$PREFIX/bin:\$PATH\"" $SHELL_RC || echo "export PATH=\"$PREFIX/bin:\$PATH\"" >> $SHELL_RC
    grep -qxF "export LD_LIBRARY_PATH=\"$PREFIX/lib:\$LD_LIBRARY_PATH\"" $SHELL_RC || echo "export LD_LIBRARY_PATH=\"$PREFIX/lib:\$LD_LIBRARY_PATH\"" >> $SHELL_RC
    grep -qxF "export PKG_CONFIG_PATH=\"$PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH\"" $SHELL_RC || echo "export PKG_CONFIG_PATH=\"$PREFIX/lib/pkgconfig:\$PKG_CONFIG_PATH\"" >> $SHELL_RC

    echo "Paths added to $SHELL_RC. Run 'source $SHELL_RC' to apply them now."
fi



