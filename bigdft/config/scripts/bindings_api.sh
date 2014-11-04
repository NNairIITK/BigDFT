#!/bin/bash
# To be called with:
#  config/scripts/bindings_api.sh src/bindings | tee src/bindings/bindings_api.h

omit=`dirname $0`/bindings_api.omit
API=`grep -rh "^ *FC_FUNC[^#]*$" $1 | cut -d'(' -f2 | cut -d',' -f1 | sort -u | grep -vf $omit`
GETS=`grep -rh "^  *GET_ATTR" $1 | cut -d'(' -f2 | sed -e "s/ //g" | cut -d',' -f1,3 --output-delimiter="_get_" | sort -u | grep -vf $omit`

echo "#ifndef BINDINGS_API_H"
echo "#define BINDINGS_API_H"
echo
echo "#undef hz"
echo
srcdir=`dirname $1`
for func in $GETS $API ; do
    func_case=`echo $func | sed -e 's/\(.\)/[\1\U\1]/g'`
    file=`grep -nr --exclude=\*interfaces.f90 --exclude=\*unused\* --exclude=\*_fake.f90 --exclude=\*_old.f90 --exclude=\*private_api.f90 --exclude=\*bindings_api.h "^ *[sS][uU][bB][rR][oO][uU][tT][iI][nN][eE] *$func_case *(" $srcdir | cut -d':' -f1,2`
    echo "/* "$func $file" */"
    python `dirname $0`/bindings_api.py $func $file
done
echo "#endif"
