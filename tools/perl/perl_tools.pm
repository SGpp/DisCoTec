#!/usr/bin/perl -w

sub max {  
 # returns the maximum of two input values 
 my $a = shift;
 my $b = shift;
 return $a > $b ? $a : $b;
}

sub min {  
 # returns the minimum of two input values 
 my $a = shift;
 my $b = shift;
 return $a < $b ? $a : $b;
}

sub isint{
  # returns true if variable is an integer number
  my $val = shift;
  return ($val =~ m/^\d+$/);
}#

sub isnumber{
  # returns true if variable is a number
  my $val = shift;
  return (($val =~ m/^((\d+)?\.)?\d+([eE][-+]\d+)?$/) or ($val=~m/^(\d+)(\.)([eE][+-]\d+)?$/) )   ;
}#

sub floor{ # round towards zero (positive numbers)
  my $val = shift;
  $val =~ /^(\d+)\.?(\d*)$/;
  return $1
}#

sub ceil{ # round towards infinity (positive numbers)
  my $val = shift;
  $val =~ /^(\d+)\.?(\d*)$/;
  return $1+1
}#

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}

sub sleep_random {
  # sleep for a short random time
  # seed for random numbers is based on 'string_in' (e.g PROBDIR)
  my $string_in = shift;
  my $len = length($string_in);
  my $sum = 0;
  for(my $i=0; $i<$len; $i++) {
    my $char = substr($string_in, $i, 1);
    $sum = $sum + int(ord($char));
  }
  #print "the sum of the ASCII-values of string is $sum\n";
  srand($sum);   #initialize random using $sum as seed
  my $sleeptime = rand(3);  #random number 0...3
  sleep($sleeptime);
}


##the following line is the end of a perl module
1;
