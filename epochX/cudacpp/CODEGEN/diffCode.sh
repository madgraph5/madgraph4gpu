#!/bin/bash
diff --no-dereference -x '*log.txt' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' $* 
