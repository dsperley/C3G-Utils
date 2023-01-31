### Loading Cairo
```
> library(Cairo)
Error: package or namespace load failed for ‘Cairo’:
 .onLoad failed in loadNamespace() for 'Cairo', details:
  call: dyn.load(file, DLLpath = DLLpath, ...)
  error: unable to load shared object '/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Cairo/libs/Cairo.so':
  dlopen(/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Cairo/libs/Cairo.so, 6): Library not loaded: /opt/X11/lib/libXrender.1.dylib
  Referenced from: /Library/Frameworks/R.framework/Versions/4.1/Resources/library/Cairo/libs/Cairo.so
  Reason: image not found
  ```
Fix: Insall XQuartz
https://stackoverflow.com/questions/38952427/include-cairo-r-on-a-mac



