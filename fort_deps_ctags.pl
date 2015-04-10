#!/usr/bin/perl
# This script analyzes fortran source code and generates a Makefile on stdout
# It requires linux ctags and takes the top level directory of the source
# tree as an optional argument
# Russell Poyner 2015

our (%prov, %uses, %file_deps);
our (%mains, @sources, @flist);
our (%make_targets, %obj_lists);
our ($fstring);

# Build the %uses hash that tells what units are referenced in each souce file
sub fortran_parser{
  my($file)=shift;
  open FORT, $file or print "couldn't open $file\n";
  my $skip = 0;
  my $iface = 0;
  my $fixedform = $file =~ /\.f$/io; # Pretty safe to assume that .f files are fixed-format
  # grep filters out un-interesting lines, speeds up parsing
  foreach my $line (grep /(interface|use|call|subroutine|function|$fstring)/i, <FORT>){
    next if($line=~/^\s*\!.*/io); # Avoid f90 comments
    next if($fixedform && $line =~ /^\S/o); # avoid fixed-format comments
    next if($line =~ /^#/o); # skip pre-processor directives
    $line =~ s/'.+'//go; # Get rid of quoted text
    $line =~ s/".+"//go;
    $line =~ s/^(.*)(!.*)/\1/o; # truncate lines after the first !
    # Add units from interface blocks to the uses list
    $iface = 0 if($line =~ /^\s*end\s+interface/io);
    if($line =~ /^\s*interface/i || $iface){
      $iface = 1;
      if($line =~ /\b(subroutine|function)\s+(\w+)/io){
        my $type = lc $1;
        my $unit = lc $2;
        push(@{$uses{$file}}, $unit);
      }
    }
    if($line =~ /\b(use|call)\s+(\w+)/io){
      $unit= lc $2;
      push(@{$uses{$file}}, $unit);
    }
    if($line=~/^(.*)\W+($fstring)(\W+.*)$/i && $1 !~ /(function|real)/io && $3 !~ /\=/o){
      push(@{$uses{$file}}, $2);
    }
  }
  close FORT;
  foreach my $key (keys %uses){
    @{$uses{$key}} = &uniquify(@{$uses{$key}});
  }
}


# Create a tree of file level dependencies for the given program unit
sub dep_tree{
  my $node = $_[0];
  $node->{'parent'} = $_[2];
  $node->{'this_file'} = $_[1];
  @{$node->{'units'}} = @{$uses{$node->{'this_file'}}};
  foreach $unit (@{$node->{'units'}}){
    next if grep /$node->{'this_file'}/, @{$prov{$unit}};
    my $pfile = @{$prov{$unit}}[0];
    unless(grep(/$pfile/, @{$node->{'files'}})) {
      push(@{$node->{'files'}}, $pfile);
      my $child = &new_node($pfile,$node->{'this_file'});
      my $ref = &dep_tree($child, $pfile, $node->{'this_file'});
      push (@{$node->{'children'}}, $ref);
    }
  }
  return $node;
}

# Print the file dependency tree
sub print_tree {
  my $node = $_[0];
  my $offset = $_[1] || 0;
  print '   ' x $offset;
  print $node->{'this_file'},"\n";
  foreach $child (@{$node->{'children'}}) {
    &print_tree($child, $offset + 1);
  }
}

# Return a reference to a new node in the dependency tree
sub new_node {
  return {
    'parent'  => $_[1]||undef,
    'this_file' => $_[0],
    'units' => [],
    'files' => [],
    'children' =>  [],
  };
}

# Return an array of unique, non-empty values
sub uniquify {
  my %saw = undef;
  grep /\S/o, grep(!$saw{$_}++, @_);
}

# Build a hash of make targets from a source dependency tree
sub targets{
  my $node = $_[0];
  my $main = $_[1];
  my $flags = '$(FC) $(FFLAGS) $(OTHERFLAGS) -c -o';
  foreach $child (@{$node->{'children'}}){
    &targets($child,$main);
  }
  return unless my $tfile = $node->{'this_file'};
  $tfile =~ s/(\.f|\.f90)$/\.o/io;
  $tfile =~ s/(.*)(\/.*)/\2/go;
  $tfile = "obj".$tfile;
  my @fdeps = map {
    s/(\.f|\.f90)$/\.o/io;
    s/(.*)(\/.*)/\2/go;
    "obj".$_;
  } grep /\w/, @{$node->{'files'}};
  my $target = "$tfile: $node->{'this_file'} @fdeps\n	$flags $tfile $node->{'this_file'}";
  $make_targets{$target} = grep /\w/o, @{$node->{'files'}};
  $make_targets{$target} += @sources unless $node->{'parent'};
  unshift @{$obj_lists{$main}}, $tfile;
}

my $base_dir = $ARGV[0] || '.';

# Use ctags to get a hash of interesting program units and the files that define them
foreach my $tagline (`ctags -R -x --languages=fortran --fortran-kinds=fspmc $base_dir`) {
  (my $unit, my $type, my $ln, my $file_name, my @rest) = split /\s+/o, $tagline;
  $mains{$unit} = $file_name if $type eq 'program';
  $unit = lc $unit;
  $type = lc $type;
  if ($prov{$unit}) {
    push @{$prov{$unit}}, $file_name;
  } else {
    $prov{$unit} = [$file_name];
  }
  push @sources, $file_name;
  push @flist, $unit if $type eq 'function';
}
@sources = &uniquify(@sources);
$fstring = join '|', &uniquify(@flist);



# Parse the source files to get dependencies
foreach my $source (@sources) {
  &fortran_parser($source);
}

# Loop over all "main" programs defined in this source tree
# and build file-based dependency tree for each
# then build target list using the trees
foreach $prog (keys(%mains)){
  my $ref = &new_node($mains{$prog});
  $ref = &dep_tree($ref, $mains{$prog}, undef);
  #print "\n==========================================================================\n\n";
  #&print_tree($ref);
  &targets($ref,$prog);
}

# Include the header
print "include make_header\n\n";

# print the compile targets
foreach $ndeps (sort {$make_targets{$a} <=> $make_targets{$b}} keys %make_targets){
  print $ndeps,"\n\n";
}

# print the object lists
foreach $prog (keys %obj_lists){
  #@{$obj_lists{$prog}} = uniquify(@{$obj_lists{$prog}});
  print 'OBJ_',$prog,' = ', join ' ', uniquify(@{$obj_lists{$prog}}),"\n\n";
}

# print the source list
print 'SRC = ',join ' ', @sources,"\n\n";

# print the linker targets
foreach $prog (keys %mains) {
  my $binfile = "bin/".$prog;
  my $objvar = '$(OBJ_'.$prog.')';
  print $binfile,': ',$objvar,"\n	",'$(LD) ',$objvar,' -o ',$binfile,' $(LDFLAGS)',"\n\n";
}

# print some additional targets
print "clean:\n	rm obj/* bin/* mod/*\n\n";
print "make:\n	\@mv Makefile Makefile.bak\n	\@./fort_deps_ctags.pl>Makefile\n\n";
print "TAGS: \$(SRC)\n	\@etags \$(SRC)\n\n";
print "tags: \$(SRC)\n	\@ctags \$(SRC)\n\n";
exit;
