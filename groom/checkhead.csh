#! /bin/csh

foreach f ( */*.c )
set strip = `echo $f | awk -F/ '{print $2}'`
head -1 $f | grep $strip > /dev/null
set h = $status
if ( $h ) then
echo $f
endif
end
