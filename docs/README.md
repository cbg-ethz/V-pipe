# Tutorials

We strongly advise our users to start discovering V-pipe by looking at the tutorials

You can find several tutorials in this directory:

## Getting V-pipe installed

- [V-pipe Installation](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_0_install.md)

## Viruses

- [V-Pipe HIV Tutorial](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_hiv.md): uses HIV test data
- [SARS-CoV-2 Tutorial](https://github.com/cbg-ethz/V-pipe/blob/master/docs/tutorial_sarscov2.md): uses SARS-CoV-2 data from a publication

## Note about the tutorials

Due to automated testing, each copy-pastable block begins with a command entering the directory and ends with one leaving the directory:

```bash
cd tutorial/work/
# do something
cd ../..
```
Of course, you don't necessarily need to do that.  You can simply remain in the working directory.

When editing files like `config.yaml`, you can use your favorite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.). By default, our tutorials use a [_heredoc_](https://en.wikipedia.org/wiki/Here_document) to make it easier to copy-paste the blocks into bash:

```bash
cat > config.yaml <<EOF
general:
    virus_base_config: 'sars-cov-2'
EOF
```
