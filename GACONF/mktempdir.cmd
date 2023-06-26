if (-e $tmpdir/$attempt.$runit) then
  @ curr = 0
  while (-e $tmpdir/$attempt.$runit.$curr)
    @ curr ++
  end
  mkdir $tmpdir/$attempt.$runit.$curr
  set dir=$tmpdir/$attempt.$runit.$curr
else
  mkdir $tmpdir/$attempt.$runit
  set dir=$tmpdir/$attempt.$runit
endif
echo $pid > $dir/proc_in_charge
