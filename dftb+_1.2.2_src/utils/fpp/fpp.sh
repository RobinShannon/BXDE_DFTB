#!/bin/sh

if [ "$1" = "-h" ]; then
    echo ""
    echo "$0 <version>"
    echo ""
    echo "Invokes afterprocessor for folding lines longer than allowed in F95"
    echo ""
    echo "<version>: general -- general folding with proper #line directives"
    echo "           nocntln -- no #line directives in continuation lines"
    echo "           noln    -- no #line directives at all."  
    echo "           noln2   -- noln + remove #line directives added by cpp"
    echo ""
    exit
fi

mypath=`dirname $0`
if [ "$1" = "noln2" ]; then
   delhash="1"
   name="noln"
else
   delhash=""
   name="$1"
fi
awk="$mypath/ff_${name}.awk"

if [ -e $awk ]; then
    if [ -z "$delhash" ]; then
      cat | awk -f $awk
    else
      cat | awk -f $awk | grep -v -e "^#"
    fi
else
    echo "Error: Awk script $awk not found" >& 2
    echo "Type '$0 -h' for help." >& 2
    exit 1
fi
