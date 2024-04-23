find . -maxdepth 1 -name "*_F_*ab1" -print | while read line; do
                 a=${line%%_*}
                 b=${a##*_}
                 c=${b//.}
                 d=${c///}
		 e=${a##*-}
                 mv $line "$e"".F.ab1";done 




find . -maxdepth 1 -name "*_R_*ab1" -print | while read line; do
                 a=${line%%_*}
                 b=${a##*_}
                 c=${b//.}
                 d=${c///}
                 e=${a##*-}
                 mv $line "$e"".R.ab1";done 
