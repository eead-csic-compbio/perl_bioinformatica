# Paul Pavlidis, 2/2000
package Stats;
require Exporter;
use POSIX;
@ISA = qw(Exporter);
@EXPORT = qw(
	     maxmin 
	     mean 
	     stdev 
	     median
	     median_old
	     quantile
	     leastsquares 
	     ranks
	     dot_product 
	     correlation_op 
	     variance 
	     tprob 
	     devprob
	     factorial 
	     lngamma 
	     binomialcoeff 
	     beta 
	     betainc  
	     betacf 
	     readlayout
	     readclass
	     factln 
	     uprob
	     fprob
	     openR);

########################
# globals
########################
my $ntop = 4;my @facttable = (1,1,2,6,24);
my @factlntable = ();
my $MINNUM = -1e100;
my $SMALL = 1.0e-30;
my $MAXNUM = 1e100;
my $MAXIT = 100;
my $EPS = 3.0e-7;
my @utable = ();

BEGIN {
  for ($i=0; $i<101; $i++) {
    $factlntable[$i] = 0;
  }
}
########################

############################
# The functions
############################

# find the min and max of an array
sub maxmin {
  my $array = shift @_;
  my ($max, $min);
  $max = $MINNUM;
  $min = $MAXNUM;
  if (scalar @$array == 0) {
    die "Empty array for Stats::maxmin\n";
  }

  foreach $value (@$array) {
    if ($value < $min) {
      $min = $value;
    } elsif ($value > $max) {
      $max = $value;
    }
  }
  return ($min, $max);
}

################################################
# Find the mean of an array of values.
sub mean {
  my $array = shift @_;
  my ($sum, $i, $value);
  if (scalar @$array == 0) {
    die "Empty array for Stats.pm::mean\n";
  }
  $sum = 0;
  $i = 0;
  foreach $value (@$array) {
    $sum+=$value;
    $i++;
  }
  return $sum / $i;
}

################################################
# Find the variance of an array of values, given the
# mean. (calculate if necessary)
sub variance {
  my ($array, $mean) = @_;
  
  # check for array of reasonable size
  if (scalar @$array < 2) {
    die "Empty or too short array for Stats.pm::stdev\n";
  }

  # calcualate mean if needed
  if (!defined $mean) {
    $mean = &mean($array);
  }

  # do stdev calc, skipping nonnumeric and missing values.
  my $i = 0;
  my $sumdev = 0;
  foreach $value (@$array) {
    $sumdev+= ($value - $mean) ** 2;
    $i++;
  }

  # recheck that array is still big enough after skipping values
  if ($i <= 1) {
    die "Empty or too short array for Stats.pm::stdev\n";
  }
  return $sumdev / ($i - 1);
}

################################################
# Find the standard deviation of an array of values, given the mean
sub stdev {
  my ($array, $mean) = @_;
  $mean = mean ($array) unless defined $mean;
  return sqrt(variance($array, $mean));
}


################################################
# Calculate the median of an array
sub median {
  my ($array) = @_;
  return quantile($array, 50);
}

sub median_old {
  my ($array) = @_;
  my @sorted;
  @sorted = sort {$a <=> $b} @$array;
  my $n = scalar @sorted;
  if ($n % 2) {
    return $sorted[floor($n/2)];
  } else {
    return ($sorted[$n/2 - 1] + $sorted[$n/2]) / 2;
  }
}

################################################ 
# Calculate any quantile of an array (quantile given as a percentage,
# rounded to an integer.). Take into account an optional array of
# weights, which must be 1 or less per datum.
sub quantile {
  my ($array, $quantile, $weights, $returndataover) = @_;
  my @sorted;
  if (@$weights) {
    if (scalar @$weights != scalar @$array) {
      die "Weights array length does not equal array length\n";
    }
    my $i;
    my @in;
    my @index;
    my $n = scalar @$array;
    for ($i=0; $i<$n; $i++) {
      push @index, [($$array[$i], $$weights[$i])];
    }

    # sort the index according to the scores
    @sorted = sort {
      $a->[0] <=> $b->[0];
    } @index;

    $quantile = (100 - $quantile)/100;

    my ($j, $totalweight);
    foreach $j (@$weights) {
      $totalweight += $j; # this is the "virtual" number of samples we have.
    }
#    print STDERR "====\n";
#    print STDERR "Total weight is $totalweight; quantile is $quantile.\n";
    if ($totalweight == 0) {
      die "Fatal: Weights totaled to zero\n";
    }
    my $k;
    my $count = 0;
    my $p = $0;
    my $quantileindex = 0;
    my $set = 0;
    foreach $k (@sorted) {
      $count+=$k->[1]; # sum the weights, which is the same as counting elements if the weights are all 1.
      $var = $count/$totalweight;
#      print STDERR "$k->[0] $k->[1] $count $var $p\n";
      if (!$set && $var > $quantile) {
	if ($p > 0) {
	    $quantileindex = $p-1; # go back one.
	  } else {
	    $quantileindex = 0;
	  }
#	print STDERR "Over: Quantile index is $quantileindex\n";
	$set++;
#	last; # comment out for debugging.
      } elsif (!$set && $var == $quantile) {
	$quantileindex = $p;
#	print STDERR "Hit: Quantile index is $quantileindex\n";
	$set++;
#	last; # comment out for debugging.
      } else {
	$p++;
      }
    }
#    print STDERR "Quantile index: $quantileindex :: Quantile value: $sorted[$quantileindex][0]\n";

    if ($returndataover) {
      my @returnvalue;
      for ($a = $quantileindex+1; $a <= $n; $a++) {
	push @returnvalue, $sorted[$sorted[$n-$a]][0];
      }
      @returnvalue = sort {$a <=> $b } @returnvalue;
      return @returnvalue;
    } else {
      return ($sorted[$quantileindex][0]);
    }
  } else { # no weights
    $quantileindex = quantileindex($array, $quantile);
    @sorted = sort {$a <=> $b} @$array;
    my $n  = scalar @$array;
    if ($returndataover) {
      my @returnvalue;
      for ($a = $n-1; $a >= $quantileindex; $a--) {
	push @returnvalue, $sorted[$a];
      }
      @returnvalue = sort {$a <=> $b } @returnvalue;
      return @returnvalue;
    } else {
      return $sorted[$quantileindex];
    }
  }
}

# return the index in an array corresponding to the quantile.
sub quantileindex {
  my ($array, $quantile) = @_;

  if ($quantile < 0 || $quantile > 100) {
    die "Quantile must be a number from 0 to 100.";
  }
  $quantile /=100; # turn into a fraction.
  my $n = scalar @$array;
  return ceil($n-1)*$quantile;
}

################################################
# Calculate the linear least squares fit for a list of x,y pairs.
sub leastsquares {
  my ($xarr, $yarr) = @_;
  my( $x, $y);
  my ($results);
  for ($i=0; $i<$#$xarr; $i++) {
    $x = $$xarr[$i];
    $y = $$yarr[$i];
    $n++;
    $sxy += $x*$y;
    $sxx += $x*$x;
    $sy += $y;
    $sx += $x;
  }
  $m = ($n*$sxy - $sy * $sx) / ($n * $sxx - $sx *$sx) ;
  $b = ($sy - $m*$sx) / $n;
  @results = ($m, $b);
  return \@results;
}

##############################################
# Dot product
sub dot_product {
  my ($array1, $array2, $dp) = @_;
  my $i, $n;
  $i = 0;
  foreach $n (@$array1) {
    $$dp += $n * $$array2[$i];
    $i++;
  }
}


################################################
# Two-sided p value (NRinC)
sub tprob {
  my ($dof, $t) = @_;
  return betainc(0.5*$dof, 0.5, $dof/($dof+$t*$t));
}

################################################
# Probability of a normal deviate. Uses an approximation that works
# for Z down to about 0.1
sub devprob {
  my ($z) = @_;
  $z = abs($z);
  my $c = $z * (1.237 + 0.0249*$z);
  return (1-sqrt(1-exp(-$c*$c)))/2;
}

################################################
# Calculate the correlation coefficient for two vectors. Pearson product moment correlation coefficient
# Assumes we can precalculate stuff about vector x.
sub correlation_op {
  my ($sumx, $sumxs, $meanx, $x, $meany, $y, $outliers) = @_; # arrays must be same size...
  my $i;
  my ($sumxy, $sumy, $sumys);
  my $numfeat = scalar @$y;
  $sumy = 0;
  $sumys = 0;
  $sumxy = 0;
  for ($i=0; $i < $numfeat; $i++) {
    $sumy += $$y[$i];
    $sumys += $$y[$i] * $$y[$i];
    $sumxy += $$y[$i] * $$x[$i];
  }
  if ($sumx*$sumx/$numfeat == $sumxs  || $sumy*$sumy/$numfeat == $sumys) {return 0;}
  return ($sumxy - ($sumx*$sumy)/$numfeat)/sqrt( ($sumxs - ($sumx * $sumx)/$numfeat) * ($sumys - ($sumy * $sumy) /$numfeat));
}



################################################
# numerical recipies in C, page 214.
sub factorial {
  my ($n) = @_;
  my $j;
  # check that $n is an integer.
  $nint = int $n;
  die "I don't think $n is an integer, you can't take the factorial of that\n" unless $nint == $n;
  
  if ($n < 0) {
    die "Attempt to calculate the factorial of a negative number\n";
  }
  if ($n > 32) {
    return exp(lngamma($n+1.0));
  }
  
  while ($ntop < $n) {
    $j=$ntop++;
    $facttable[$ntop] = $facttable[$j]* $ntop;
  }
  return $facttable[$n];
}

################################################
# numerical recipes in C, page 214
sub lngamma {
  my ($n) = @_;
  @cof = (76.18009172947146,
	  -86.50532032941677,
	  24.01409824083091,
	  -1.231739572450155,
	  0.120865097386617e-2,
	  -0.5395239384953e-5);
  $y = $n;
  $x = $n;
  $temp = $x + 5.5;
  $temp -= ($x + 0.5)*log($temp);
  $ser = 1.000000000190015;
  for ($j=0; $j<=5; $j++) {
    $ser+= $cof[$j]/++$y;
  }
  return - $temp + log(2.5066282746310005*$ser/$x);
}

################################################
# numerical recipes in C, page 215
# returns ln ($n!)
sub factln {
  my ($n) = @_;
  die "Attempt to calculate the factorial of a negative number\n" if ($n < 0);
  return 0 if ($n <= 1);
  if ($n <= 100) {
    if (!$factlntable[$n]) {
      $factlntable[$n] = lngamma($n+1);
    }
    return $factlntable[$n];
  } else {
    return lngamma($n+1);
  }
}

################################################
# numerical recipes in C, page 215
sub binomialcoeff {
  my ($n, $k) = @_;
  return floor(0.5 + exp(factln($n) - factln($k) - factln($n-$k)));
}

################################################
# numerical recipes in C, page 216
sub beta {
  my ($z, $w) = @_;
  return exp (lngamma($z) + lngamma($w) - lngamma($z+$w));
}

################################################
# numerical recipes in C, page 227
sub betacf {
  my ($a, $b, $x) = @_;
  my ($m, $m2);
  my ($aa, $c, $d, $del, $h, $qab, $qam, $qap);
  $qab = $a + $b;
  $qap = $a + 1.0;
  $qam = $a - 1.0;
  $c = 1.0;
  $d = 1.0 - $qab*$x/$qap;
  if (abs($d) < $SMALL) { $d = $SMALL; }
  $d = 1.0/$d;
  $h = $d;
  for ($m=1; $m<=$MAXIT; $m++) {
    $m2=2*$m;
    $aa=$m*($b-$m)*$x/(($qam+$m2)*($a+$m2));
    $d=1.0+$aa*$d;
    if (abs($d) < $SMALL) { $d = $SMALL; }
    $c = 1.0 + $aa/$c;
    if (abs($c) < $SMALL) { $c = $SMALL; }
    $d = 1.0/$d;
    $h *= $d * $c;
    $aa = -($a + $m)*($qab+ $m)*$x/(($a + $m2)*($qap+$m2));
    $d=1.0+$aa*$d;
    if (abs($d) < $SMALL) { $d = $SMALL; }
    $c=1.0+$aa/$c;
    if (abs($c) < $SMALL) { $c = $SMALL; }
    $d=1.0/$d;
    $del = $d*$c;
    $h *= $del;
    if (abs($del - 1.0) < $EPS) { last; }
  }
  if ($m > $MAXIT) { die "Could not calculate beta continued fraction value\n"; }
  return $h;
}

################################################
# incomplete beta function
sub betainc {
  my ($a, $b, $x) = @_;
  my ($bt);
  if ($x < 0 || $x > 1) { die "Illegal value $x for x in incomplete beta function\n";  }
  if ($x == 0 || $x == 1.0) {
    $bt  = 0.0; 
  }  else {
    $bt = exp (lngamma($a+$b) - lngamma($a) - lngamma($b) + $a*log($x) + $b*log(1-$x));
  }
  if ($x < ($a + 1)/($a + $b + 2)) {
    return $bt*betacf($a, $b, $x)/$a;
  }  else {
    return 1.0 - $bt * betacf($b, $a, 1-$x)/$b;
  }
}

################################################
# nrc pg 229.
sub fprob {
  my ($f, $v1, $v2)  = @_;
  return betainc($v2/2, $v1/2, $v2/($v2 + $v1*$f));
}

############ return an array of ranks corresponding to an array of data.
sub ranks {
  my ($data, $result)  = @_;
  my @datasort;
  my %dataloc;
  my $index = 0;
  my $m;

  # build a hash telling where each value is contained in the original data set.
  foreach $m (@$data) {
    push @{$dataloc{$m}} , $index; # need to use array because of ties.
    $index++;
  }
  
  @datasort = sort { $a <=> $b } @$data;
  my $rank = 1; # be careful; if ranks start at  0, the sum of the ranks will be n(n-1)/2 instead of n(n+1)/2
  foreach $m (@datasort) {
    if (scalar @{$dataloc{$m}} > 1) { # tie
      $spreadrank = 0;
      foreach $occurrence (@{$dataloc{$m}}) {
	$spreadrank += $rank;
	$rank++;
      }
      $spreadrank = $spreadrank / scalar @{$dataloc{$m}};
      $count = 0;
      foreach $occurrence (@{$dataloc{$m}}) {
	$$result[$occurrence] = $spreadrank; # same value, same rank.
	if ($count > 0) {
	  shift @datasort;
	}
	$count++;
      }
    } else {
      $$result[$dataloc{$m}->[0]] = $rank; 
      $rank++;
    }
  } 
}

# U distribution probability
sub uprob  {
  my ($u, $n1, $n2) = @_;
#  if ( ( $n1 < 20 && $n2 < 40 ) || ( $n1 < 40 && $n2 < 20)) {
    # lookup
    
#  } else {
    # use normal approximation. Not accurate for small ns.
    my $mean = $n1 * $n2 / 2;
    my $stderr = sqrt($n1*$n2*($n1+$n2 + 1)/12);
    my $z = abs ($u - $mean) / $stderr;
    return tprob(10000, $z); # tprob for infinite degrees of freedom is what we're supposed to use.
#  }
}


###############################################
# REad a layout file
###############################################
sub readlayout {
  my ($layout, $cat) = @_;
  open (IN, "<$layout") or die "Could not open layout $layout\n";
  my ($m, $k, $q, $category, @dat, $levelname);
  $/= "=";
  <IN>;
  my %used;
  my @categories;
  my %levelused;
  while (<IN>) {
    chomp;
    s/\cM//g;
    ($category, @dat) = split "\n", $_;
    $category =~ s/=//;
    print STDERR "Variable: $category";
    push @categories, $category;
    %used = {};
    
    foreach $m (@dat) {
      next if $m eq ""; # allow blank lines.
      if($m =~ /\%(.+)/) {
	$levelname = $1;
	print STDERR "\n\tLevel: $levelname\n";
	%levelused = {};
      } elsif ($m =~ /^[0-9]+$/) {
	if ($used{$m}) {
	  die "The same index $m was used twice in this variable ($category)\n";
	}
	if ($levelused{$m}) {
	  die "The same index $m was used twice in this level ($levelname for $category)\n";
	}
	push @{$cat->{$category}->{$levelname}}, $m;
	$used{$m}++;
	$levelused{$m}++;
	print STDERR "\t$m ";
      } else {
	die "Illegal line in the layout file\n";
      }
    }
    print STDERR "\n";
  }
  close IN;
  $/="\n"; # very important.
  $numcategories = scalar keys %levelused;

  # error checking
  my %replicate;
  my %numwithcats;
  die "No variables found in the layout file\n" if (scalar %$cat == 0);
  foreach $category (keys %$cat) {
    die "No levels found for $category\n" unless (scalar keys %{$cat->{$category}} > 0);
    foreach $level (keys %{$cat->{$category}}) {
      die "No data found for $level::$category\n" unless scalar (@{$cat->{$category}->{$level}});
      foreach $q (@{$cat->{$category}->{$level}}) {
	push @{$numwithcats{$q}}, $level; # associate categories with this trial.
      }
    }
  }

  # check
  if ($debug) {
    foreach $k (sort keys %numwithcats) {
      $catlist = join " ", @{$numwithcats{$k}};
      push @{$replicate{$catlist}}, $k; # reverse the hash...
    }
    foreach $k (sort keys %replicate) {
      foreach $m (@{$replicate{$k}}) {
	print STDERR "$m $k\n";
      }
    } 
  } # debug code

  return $numcategories;
} # readlayout

#######################################################################
# readclass: read a class file, return same data structure as readlayout
sub readclass {
  my ($classfile, $cat) = @_;
  open (IN, "<$classfile") or die "Could not open $classfile\n";
  ($category = $classfile) =~ s/\.txt$//;
  $category =~ s/(.+)\///;
  print STDERR "$category\n";
  <IN>;
  my $m = 0;
  my %numcat;
  while ($line = <IN>) {
    chomp $line;
    $line =~ s/\cM//;
    if ($line eq "") {
      next; # skip blank lines
    } else {
      ($sample, $class) = split "\t", $line;
      push @{$cat->{$category}->{$class}}, $m;
      $numcat{$class}++;
      $m++;
    }
  }
  foreach $m (keys %{$cat}) {
    foreach $k (keys %{$cat->{$m}}) {
      foreach $p (@{$cat->{$m}->{$k}}) {
	print STDERR "$p\t$k\n";
      }
    }
  }
  close IN;
  $numcategories = scalar keys %numcat;
  return $numcategories;
}  # readclass

1;

# open R on a pipe from stdin.
sub openR {
  my ($handle) = @_;
  open ($handle, "|R CMD BATCH --no-save --no-restore --slave") or die "Could not open R: $!\n";
}
