# V-pipe
Vpipe Jekyll custom site draft version

## Navigation
Navigation items can be added/removed in the nav_items.yml file in the _data folder. The display logic is in navigation.html located in the _includes folder.
Note: [Footer](#footer) is independent.

## Hero carousel
Content can be added/edited in the heroCarousel collection. 
* Ideally, images should be sized at 1460x872 (in a 5:3 aspect ratio) to optimize for retina displays, with a minimum size of 730x438.

## Get started with V-pipe
Content can be added/edited in the getStartedColumns collection

## V-pipe in numbers
This is currently hard-coded in vpipe_number.html in _includes

## Why choose V-pipe
Content can be added/edited in the whyVpipe collection

## Meet the V-pipe team
**To add an event** add/remove an entry in the [events.yml](https://github.com/cbg-ethz/V-pipe/blob/gh-pages/_data/events.yml) file in the _data folder. 
The logic will check the data against the current date and put them either in the upcoming- or past events.

## V-pipe funding
Content is hard-coded in the funding.html file in the _includes folder

## Adding articles to Literature

**To add a publication** to the "[Literature](https://cbg-ethz.github.io/V-pipe/literature/)" page:
1) get your BIB.  Most pre-print and pub servers have a "Download BIB" button. (otherwise for details from PubMed, see `_literature/HOWTO.md`).
2) add the bib entry to correct file in `_literature/<<section>>.bib` depending on your contribution to V-Pipe:
   i.e. `<<section>>.bib` as `primary.bib`, `general.bib` or `components.bib`.
3) navigate to `_literature/`
4) run `make` to build all the new resources, commit all.
   > [!NOTE]  
   > Note you may need to install "pandoc" and "bibutils" before.

## Footer
Content can be added/removed in the _footer collection and the footer.html in the _includes folder.
Note: [Navigation](#navigation) is indepedent
