# Tools for working with sequencing data
### used for 
+ creating panel designs
+ add transcriptional and mutational data

### available as code best used from jupyter notebooks or docker container
+ #### for jupyter:
   * clone the repository into \<your folder\> and move into pyseq:
      + `cd <your_folder> && git clone git@github.com:Mar111tiN/pyseq.git && cd pyseq`

+ #### for docker:
   * see on [Docker Hub](https://hub.docker.com/repository/docker/martin37szyska/primertools)
   * `docker run -p <local-IP>:8888 -v $(pwd):/home/martin/work -v <your_static_folder>:/home/martin/static martin37szyska/primertools:<TAG>`
   * store your own code in `<mounted_volume>/code` and import directly (folder code automatically added to PYTHONPATH)
   * example notebooks are available in nb_templates
