<?php

echo 'Hello world';
$myvar=shell_exec("cat numeration.txt");
//echo $myvar;
$result=shell_exec("./script01.sh ".$myvar);
echo $result

?>
