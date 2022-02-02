#!/bin/bash
diff --no-dereference -x '*log.txt' -x 'nsight_logs' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' -x 'build.*' -x '.build.*' -x '*~' -x 'include' $* 
