#!/bin/bash

var=$(($(($(ls | wc | awk '{print $1}')-107))/3))
echo $var
