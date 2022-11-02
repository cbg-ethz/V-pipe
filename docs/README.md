# Tutorials

You can find two tutorials in this directory:

- [tutorial_hiv.md](tutorial_hiv.md): uses HIV test data
- [tutorial_sarscov2.md](tutorial_sarscov2.md): uses SARS-CoV-2 data from a publication

## Note about the tutorials

Due to automated texting, each copy-pastable block begins with a command entering the directory and ends with on leaving the directory:

```bash
cd tutorial/work/
# do something
cd ../..
```
Of course you don't necessarily need to do that.  You can simply remain in the working directory.

When editing files like `config.yaml`, you can use your favorite editor (`vim`, `emacs`, `nano`, [butterflies](https://xkcd.com/378/), etc.). By default our tutorials use a [_heredoc_](https://en.wikipedia.org/wiki/Here_document) to make it easier to copy-paste the blocks into bash:

```bash
cat > config.yaml <<EOF
general:
    virus_base_config: 'sars-cov-2'
EOF
```
