#!/usr/local/bin/perl 
$cpp = '/usr/lib/cpp  ';
$cc = 'mpicc ';

while ( $_ = shift ){
    s/>/\\>/g;
    s/</\\</g;
    if(/^-D/) { $cpp=$cpp.' '.$_ ; next ;}
    if(/^-I/) { $cpp=$cpp.' '.$_ ; next ;}
    if(/\.c$/) { $cpp=$cpp.' '.$_ ; ($file = $_)=~s/\.c$/.i/ ; next ;}
    $cc = $cc.' '.$_ ;
}
$cc = $cc.' '.$file;
$cpp = $cpp.' '.$file;

print $cpp."\n" ;
system($cpp);
print $cc."\n" ;
system($cc);
system("rm ".$file);
