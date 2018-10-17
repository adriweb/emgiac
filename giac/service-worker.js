/* global self, caches*/

/**
* Modified from http://dev.geogebra.org/trac/browser/trunk/geogebra/web/dev/src/main/resources/org/geogebra/web/worker_template.js?rev=50858 Designed by Gabor, made by GWT, thanks to Zbyněk Konečný 
*/
var urlsToCache = {"urls_to_cache":["giac.js","giacwasm.js","giacwasm.wasm","xcas.js","xcasfr.html","xcasen.html","bootstrap/bootstrap.min.css","bootstrap/bootstrap-theme.min.css","fonts/glyphicons-halflings-regular.woff2","codemirror/codemirror.css","codemirror/codemirror.js","codemirror/dialog.css","codemirror/dialog.js","xcasmode.js","python.js","codemirror/matchbrackets.js","FileSaver.js","w3data.js","menufr.js","menuen.js","logo.png","config.png","longhelp.js","longhelp_en.js","giacworker.js","codemirror/show-hint.js","codemirror/show-hint.css","codemirror/search.js","codemirror/searchcursor.js","codemirror/jump-to-line.js","algoseconde.html","mathjax/MathJax.js?config=TeX-AMS_CHTML"],
		   "unique_id":"#1.4.9-61"};

self.addEventListener('install', function(event) {
    "use strict";
    console.log("install");
    event.waitUntil(
        caches
            .open(urlsToCache.unique_id)
            .then(function(cache) {
                console.log('[install] Caches opened, adding Xcas for Firefox files to cache');
                return cache.addAll(urlsToCache.urls_to_cache);
            })
            .then(function() {
                console.log('[install] All required resources have been cached');
            })
    );
});

self.addEventListener('fetch', function(event) {
    "use strict";
    console.log('fetch');
    event.respondWith(
        caches.match(event.request)
            .then(function(response) {
                if (response) {
                    return response;
                }
                return fetch(event.request);
            }).catch(function(reason) {
                    console.log(reason);
            })
    );
});

self.addEventListener('activate', function(event) {
    "use strict";

    console.log("activate");

    event.waitUntil(
        caches.keys()
            .then(function(cacheNames) {
                return Promise.all(
                    cacheNames.map(function(cacheName) {
                        if (urlsToCache.unique_id !== cacheName) {
                            console.log("deleting from cache " + cacheName);
                            return caches.delete(cacheName);
                        }
                    })
                );
            })
    );
});
