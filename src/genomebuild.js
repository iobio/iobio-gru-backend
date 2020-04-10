const Router = require('koa-router');
const { dataPath } = require('./utils.js');

var sqlite3 = require('sqlite3').verbose();
var db = new sqlite3.Database(dataPath('genomebuild/genomebuild.db'));

const router = new Router();

router.get('/', async (ctx) => {  

  // Get all species
  var speciesSql = "SELECT * from species";
  var species = [];

  return new Promise((resolve, reject) => {
    db.all(speciesSql,function(err,speciesRows){ 

      if (err) reject(err);

      if (speciesRows && speciesRows.length > 0) {
        speciesRows.forEach(function(speciesRow) {
          species.push(speciesRow);
        });

        // Get all genome builds and map them back to parent species
        var buildSql = "SELECT * from genomebuild ";
        var genomeBuilds = [];
        db.all(buildSql,function(err,buildRows){ 

          if (err) reject(err);

          var genomeBuild = {};
          var buildMap = {};
          if (buildRows != null && buildRows.length > 0) {
            for (var i = 0; i < buildRows.length; i++) {
              genomeBuild = buildRows[i];    
              genomeBuilds.push(genomeBuild);

              var builds = buildMap[genomeBuild.idSpecies];
              if (builds == null) {
                builds = [];
                buildMap[genomeBuild.idSpecies] = builds;
              }
              builds.push(genomeBuild);
            }
            species.forEach(function(species) {
              species['genomeBuilds'] = buildMap[species.id];
            });

            // Get all references and map them back to parent genome build
            var referenceMap = {};
            var refSqlString = "SELECT * from reference";
            db.all(refSqlString,function(err,refRows){ 

              if (err) reject(err);

              if (refRows != null && refRows.length > 0) {
                refRows.forEach(function(refRow) {
                  var references = referenceMap[refRow.idGenomeBuild];
                  if (references == null) {
                    references = [];
                    referenceMap[refRow.idGenomeBuild] = references;
                  }
                  references.push(refRow);
                });
                genomeBuilds.forEach(function(genomeBuild) {
                  genomeBuild['references'] = referenceMap[genomeBuild.id];
                });

                // Get all genome build resources and map them back to parent genome build
                var resourceMap = {};
                var resourceSqlString = "SELECT * from genomeBuildResource";
                db.all(resourceSqlString,function(err,resourceRows){ 

                  if (err) reject(err);

                  if (resourceRows != null && resourceRows.length > 0) {
                    resourceRows.forEach(function(resourceRow) {
                      var resourcees = resourceMap[resourceRow.idGenomeBuild];
                      if (resourcees == null) {
                        resourcees = [];
                        resourceMap[resourceRow.idGenomeBuild] = resourcees;
                      }
                      resourcees.push(resourceRow);
                    });
                    genomeBuilds.forEach(function(genomeBuild) {
                      genomeBuild['resources'] = resourceMap[genomeBuild.id];
                    });
                  }

                // Get all genome build aliases and map them back to parent genome build
                var aliasMap = {};
                var aliasSqlString = "SELECT * from genomeBuildAlias";
                  db.all(aliasSqlString,function(err,aliasRows){ 

                    if (err) reject(err);

                    if (aliasRows != null && aliasRows.length > 0) {
                      aliasRows.forEach(function(aliasRow) {
                        var aliases = aliasMap[aliasRow.idGenomeBuild];
                        if (aliases == null) {
                          aliases = [];
                          aliasMap[aliasRow.idGenomeBuild] = aliases;
                        }
                        aliases.push(aliasRow);
                      });
                      genomeBuilds.forEach(function(genomeBuild) {
                        genomeBuild['aliases'] = aliasMap[genomeBuild.id];
                      });

                      ctx.set('Content-Type', 'application/json');
                      ctx.set('Charset', 'utf-8')
                      ctx.body = ctx.query.callback + '(' + JSON.stringify(species) +');';
                      resolve();
                    }
                  });


                });

              }
            });

          } 
        });      
      }
    });
  });
});

module.exports = router;
