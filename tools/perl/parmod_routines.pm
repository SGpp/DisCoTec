#!/usr/bin/perl -w


sub set_scanscript_lib{
  # to be called in GENE directory (by newprob):

  # copies the scanscript from /tools/perl to /tools
  # and modifies the copy to include local gene tools path
  # (useful for scanscript --mks use in output directory, for example)
  my $cwd = cwd();
  my $scanscriptdir="$cwd/tools/perl";
  my $sfile="$cwd/tools/scanscript";
  copy("$scanscriptdir/scanscript","$sfile");
  my $line='';
  my $newline='';
  open(SFILE, "$sfile");
  my @content=<SFILE>;
  close(SFILE);

  open(SFILE, ">$sfile");
  foreach $line(@content) {
     if ($line =~/use\s*lib.*\/tools\/perl/){
        $newline="use lib \"$scanscriptdir\";\n";
        if ($line !~/$newline/){
           $line=$newline;
           print "set lib path in ./tools/scanscript\n";
        }
     }
     print SFILE $line;
  }
  close(SFILE);
}

#########################################################################################
### subroutines that can also be used in the testsuite:
#########################################################################################

sub read_entry {
  # returns the value of parameter $entry in input file
  my $file = shift;
  my $entry = shift;
  my $line;
  open(FH,"<$file");
  while ($line=<FH>) {
     if ($line =~ /^\s*$entry\s*=\s*(.*)/) {
        close(FH);
        return($1);
     }
  }
  close(FH);
  return("");
}

sub read_entry_in_nmllist {
  # returns the value of parameter $entry in the namelist of the input file
  my $file = shift;
  my $namelist = shift;
  my $entry = shift;
  my $line;
  my $in_namelist = 0;
  open(FH,"<$file");
  while ($line=<FH>) {
     if ($line =~ /^\s*$namelist\s*/){
        #found the (uncommented) namelist and starts to look for entry
        $in_namelist = 1;
     }elsif ($in_namelist == 1) {
        #check for entry in current line:
        if ($line =~ /^\s*$entry\s*=\s*(.*)/) {
            close(FH);
            return($1);   
        }elsif ($line =~ /^\/\s*(.*)/) {
           #found the '/' symbol, i.e. finished
           #the namelist
           $in_namelist = 0;
        }
     }
  }
  close(FH);
  return("");
}

sub add_entry_if_missing{
  #  do nothing, if parameters exists in input file,
  #  add the parameter $entry to $namelist otherwise
  #  e.g. add_entry_if_missing("./parameters","&general","calc_dt",".true.");  
  my $filename = shift;
  my $namelist = shift;
  my $entry = shift;
  my $newvalue = shift;
  my $line;
  if (read_entry("$filename","$entry") eq ""){
     open(DAT,"$filename");
     my @lines=<DAT>;
     close(DAT);
     open(DAT,">$filename");
     foreach $line(@lines){
        if ($line =~ /^\s*$namelist\s*/){
           $line="$line"."$entry = $newvalue\n";
        }
        print DAT $line;
     }
     close DAT ;
     return("a"); #added
  }else{
    return("n"); #nothing done
  }
  print "error in set_entry : $entry";
  return(0);
}



sub comment_entry {
  # comment out line containing keyword
  my $filename = shift;
  my $keyword = shift;
  my $line;
  open(DAT,"$filename");
  my @entry=<DAT>;
  close(DAT);
  open(DAT,">$filename");
  foreach $line(@entry){
     if ($line =~ /^\s*$keyword\s*(.*)/) {
     # only uncommented lines count here 
        $line="!$line";
     }
     print DAT $line;
  }
  close(DAT);
}

sub uncomment_entry {
  # uncomment line containing keyword
  # NOTE: all "!" at the beginning of the line are erased
  my $filename = shift;
  my $keyword = shift;
  my $verbose= 0;
  $verbose=shift;
  my $line;
  open(DAT,"$filename");
  my @entry=<DAT>;
  close(DAT);
  open(DAT,">$filename");
  foreach $line(@entry){
     if ($line =~ /^(\s*!\s*)*($keyword\s*(.*))/) {
        #old: removes first character
        #$line=substr($line,1);

        #new: removes all ! in the beginning of a line
        if ($verbose==1){
           print "$filename\n";
           print "uncommented line $2\n";
        }
        $line = "$2\n";
     }
     print DAT $line;
  }
  close(DAT);
}


sub delete_entry{
  # delete the line containing keyword
  my $filename = shift;
  my $keyword = shift;
  my $line;
  open(DAT,"$filename");
  my @entry=<DAT>;
  close(DAT);
  open(DAT,">$filename");
  foreach $line(@entry){
     if ($line =~ /^\s*$keyword\s*(.*)/) {
        $line="";
     }   
     print DAT $line;
  }
  close(DAT);
}

sub set_entry{
  # replace an entry in specified namelist if it exists, 
  # add to input namelist otherwise
  my $filename = shift;
  my $namelist = shift;
  my $entry = shift;
  my $newvalue = shift;
  my $line;

  if (read_entry_in_nmllist("$filename","$namelist","$entry") eq ""){
     open(DAT,"$filename");
     my @lines=<DAT>;
     close(DAT);
     open(DAT,">$filename");
     foreach $line(@lines){
        if ($line =~ /\s*$namelist\s*/){
           $line="$line"."$entry = $newvalue\n";
        }
        print DAT $line;
     }
     close DAT ;
     return("a"); #added
  }else{
     replace_entry("$filename","$entry", "$newvalue");
     return("r"); #replaced
  }
  print "error in set_entry : $entry";
  return(0);
}


sub set_entry_all_files_in_dir{
# set an entry in all parameter files in the input directory
  my $dir = shift;
  my $namelist = shift; 
  my $par_name = shift;
  my $new_value = shift;
  my $silent = shift;
  my $elem = "";
  if ("$silent" ne "silent"){
     print "set parameter \"$par_name\" in namelist \"$namelist\" to \"$new_value\"\n";
     print "in all parameter files found in directory:\n $dir\n";
  }
  opendir(PATH,"$dir");         
  my @entry=readdir(PATH);
  closedir(PATH);
  foreach $elem(@entry){
     if($elem=~/parameters/){
        my $parfile = "$dir"."/$elem";
        set_entry("$parfile","$namelist","$par_name","$new_value")  ;
     }
  }
}


sub replace_entry {
  # replace entry or commented entry in file
  my $filename = shift;
  my $keyword = shift;
  my $newentry = shift;
  my $line;
  open(DAT,"$filename");
  my @entry=<DAT>;
  close(DAT);
  open(DAT,">$filename");
  foreach $line(@entry){
     if ($line =~ /^!*\s*$keyword\s*=\s*(.*)/) {
        $line="$keyword = $newentry\n";
     }
     print DAT $line;
  }
  close(DAT);
}


#debug routine...
sub print_scan_namelist{
  $f = shift;
  $npps = read_entry("$f","n_procs_sim");
  $nps  = read_entry("$f","n_parallel_sims");
  $pid  = read_entry("$f","par_in_dir");
  $sd  =  read_entry("$f","scandims");
  print "n_procs_sim = $npps\n";
  print "n_parallel_sims = $nps\n";
  print "scan_dims       = $sd\n";
  print "par_in_dir      = $pid\n";
}

#####################################################################################
################################ special scanscript routines ########################
#####################################################################################


sub ReadValue{
# read an entry of the input parameter file 
# arguments are:
  my $dat = shift;     #file name
  my $entry = shift;   #name of the paramter to be read, keyword
# $i   : read the $i'th occurence of the parameter $entry  
# the $all switch allows to read the whole line or skip comments like "!scan: XX"  
  my $i = 1;
  my $all = 0;          
  $i=shift;            #index of keyword to be read (i'th occurence will be read)
  $all = shift;        #if ($all=1) return everything after "=", (all=0): skip comments

  open(FILE,"<$dat");
  while(<FILE>){
    if ($all==0){   #read without comment
      if ($_ =~ /^\s*$entry\s*=\s*(\S+)\s*\!?/){
         if ($i==1){
            close(FILE);
            return($1);
         }else{
            $i--;
         }
       }
    }else{ #read all
       if ($_ =~ /^\s*$entry\s*=(.*)/) {
         if ($i==1){
            close(FILE);
            return($1);
         }else{
            $i--;
         }
       }  
   }
 }
 close(FILE);
 #print "$entry not found by routine ReadValue in file \n $dat \n";
 return(0);
}
####################################################################################


sub ChEntry{
#   change value of parameter entry in the input parameter file
#   conserves comments like "!scan: 0.1,0.1,1" at the end of the line
    my $file = shift;
    my $entry = shift;
    my $ch = shift;
    my $spec=shift;  # the $spec'th ocurrence of $entry will be changed
    my $silent=shift; # no output in silent mode
    my $i = $spec;
    my $line;
    open(FILE,"$file")|| die "couldn't open dir.";
    my @entryarr=<FILE>;
    close(FILE);    
    my $complex=0;
    if($entry=~/(.*)_scan_(re|im)/){
        $entry=$1;
        if($2=~/re/){
           $complex=1;
        }else{
           $complex=2;
        }
    }
    foreach $line(@entryarr){
        if ($line =~ /^\s*$entry\s*=.+/) {
           if($i==1){
              #if ($line =~ /!scan(\w*):((re|im)?:){1}/) {
              #   print"$entry $spec = $ch\n";
              #   $line=Replace($line,$entry,$ch,$complex);
              #}else{
               #$ch=$ch."\n"; #newline ist fuer strings wie diagdir und logicals .t. wichtig
               if ("$silent" ne "silent"){
                  print"   $entry $spec = $ch\n";
               }
               $line=Replace($line,$entry,$ch,$complex);
           }
           $i--;   #only edit $spec th occurence of $entry
        }   
    }
    open(FILE,">$file");
    print FILE @entryarr;
    close(FILE);    
    if (($i>0)&&("$silent" ne "silent")) {
        print"parameter $entry($spec) in $file not found in ChEntry \n";
    }
    return();
}
####################################################################################



sub Replace{
#   replaces a parameter value, used in the ChEntry routine  
  my $line=shift;
  my $entry=shift;
  my $ch=shift;
  my $complex=shift;
  #   general routine, can handle 
  #   "$ch !scan...\n" and "$ch\n"
  #   therefore it is important to remove \n from $line and $ch and add \n at the end
  chomp($line);
  chomp($ch);
  if($complex==0){
     $line=~ s/\s*$entry\s*=\s*\S+\s*/$entry = $ch /; 
  }elsif($complex==1){ #real part
     $ch=~/\s*(\S*)\s*/;
     $ch=$1;
     $line=~ s/\s*$entry\s*=\s*\(\s*\S+\s*,\s*(\S+)\s*\)\s+/$entry = ($ch,$1) /;
  }else{            #immaginary part
     $ch=~/\s*(\S*)\s*/;
     $ch=$1;
     $line=~ s/\s*$entry\s*=\s*\(\s*(\S+)\s*,\s*\S+\s*\)\s+/$entry = ($1,$ch) /;
  }
  #   add \n at the end
  $line = $line."\n";
  return($line);
}
####################################################################################


##the following line is the end of a perl module
1;
