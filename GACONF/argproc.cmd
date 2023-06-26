if ($#argv>0) then
  foreach a ($argv)
    if ( !($a == $a:s/=//) ) then
       eval set $a
    endif
  end
endif
