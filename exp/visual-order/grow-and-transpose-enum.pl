#!/usr/bin/perl

use strict;
use warnings;

my @matrix;
my $max_w = 0;

while (<>) {
  s/\[.*?\]//g;
  my @nums = /(\d+)/g;
  next unless @nums;

  my $code_str = join("0", map { "1" x $_ } @nums);
  my @bits = split(//, $code_str);

  push @matrix, \@bits;
  $max_w = @bits if @bits > $max_w;
}

for my $i (0 .. $max_w - 1) {
  print join(" ", map { $_->[$i] // " " } @matrix) . "\n";
}
