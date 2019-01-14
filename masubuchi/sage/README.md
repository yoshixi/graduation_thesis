
docker command

```
$ docker build -t sage .
$ docker run --rm -it -v $(pwd):/usr/src/app sage bash

# sageはcontainerから直接install
# installにめっちゃ時間かかる()
# $ sudo apt-get install sagemath
```

[sage tutorial](http://doc.sagemath.org/html/ja/tutorial/index.html)
