#!/bin/sh


ipcs -s | awk ' $3 == "'$LOGNAME'" {print $2, $3}' | awk '{ print $1}' | while read i; do ipcrm -s $i; done

ipcs -s


rm -rf /tmp/socket*

./cluster
