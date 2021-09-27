#!/bin/bash
diff -x '*log.txt' -x '*.o' -x '*.o.*' -x '*.a' -x '*.exe' -x 'lib' $* 
