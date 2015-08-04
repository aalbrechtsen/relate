.onLoad=function(libname, pkgname)
{
        library.dynam("Relate", pkgname, libname)
}

.onUnload=function(libpath)
{
        library.dynam.unload("Relate", libpath)
}
