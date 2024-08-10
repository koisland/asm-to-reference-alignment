### Changed
* Fixed alignment rule logs not being directed to stderr.
* One assembly per sample.
* No more `tbl`. Just input sample and asm to remove `pandas` dependency.
```yaml
sm:
  HG002: .test/HG002_1.fa.gz
  HG00733: .test/HG00733_1.fa.gz
```
* Removed everything but asm-ref-alignment and Saffire
* Changed benchmarks to output in benchmarks directory.
* Added option to specify directories.
* Remove secondary alignment.
