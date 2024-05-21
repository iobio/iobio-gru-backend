var async = require('async');
const Router = require('koa-router');
const router = new Router();

let _db;
function getDb() {
  if (!_db) {
    const sqlite3 = require('sqlite3').verbose();
    const { dataPath } = require('./utils.js');
    _db = new sqlite3.Database(dataPath('geneinfo/gene.iobio.db'));
  }
  return _db;
}

getGenesInClause = function(genes) {
  sqlString = " in (";
  let firstTime = true;
  genes.forEach(function(gene) {
    if (firstTime) {
      firstTime = false;
    } else {
      sqlString += ",";
    }
    sqlString += "\""+gene+"\" ";
  }) 
  sqlString += ")";
  return sqlString
}

router.get('/api/gene/:gene', async (ctx) => {  
  var source = ctx.query.source;
  var species = ctx.query.species;
  var build = ctx.query.build;
  if (source == null || source == '') {
    source = 'gencode';
  } 
  var geneSqlString = "SELECT g.*, ";
  geneSqlString += "gs.gene_symbol, ";
  geneSqlString += "GROUP_CONCAT(distinct ga.alias_symbol) AS aliases ";
  geneSqlString += "FROM genes g ";
  geneSqlString += "LEFT OUTER JOIN gene_symbol gs on gs.gene_symbol = g.gene_symbol ";
  geneSqlString += "LEFT OUTER JOIN gene_alias  ga on ga.gene_symbol = g.gene_symbol and ga.alias_symbol != g.gene_name ";
  geneSqlString += "WHERE gene_name like \""+ctx.params.gene+"\" ";
  geneSqlString += " AND source = \""+source+"\"";
  if (species != null && species != "") {
    geneSqlString  += " AND species = \""+species+"\"";
  }
  if (build != null && build != "") {
    geneSqlString  += " AND build = \""+build+"\"";
  }

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(geneSqlString,function(err,rows){ 
      var gene_data = {};
      var transcript_ids = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          gene_data = rows[i];    
          if (gene_data.hasOwnProperty("transcripts") && gene_data.transcripts != null && gene_data.transcripts != "") {
            transcript_ids = transcript_ids.concat(JSON.parse(gene_data['transcripts']));
          }       
        }
      } 
          
      async.map(transcript_ids,      
        function(id, done){      
          var source = ctx.query.source; 
          if (source == null || source == '') {
            source = 'gencode';
          }
          var sqlString = "";

          sqlString =  "SELECT t.* ";
          sqlString += "from transcripts t ";
          sqlString += "WHERE t.transcript_id=\""+id+"\" "
          sqlString += " AND t.source = \""+source+"\"";
          if (species != null && species != "") {
            sqlString  += " AND t.species = \""+species+"\"";
          }
          if (build != null && build != "") {
            sqlString  += " AND t.build = \""+build+"\"";
          }        
          db.all(sqlString,function(err,rows){    

            if (err) reject(err);

            if (rows != null && rows.length > 0) {
              rows[0]['features'] = JSON.parse(rows[0]['features']);
            } else {
              rows[0]['features'] = [];
            }   
            done(null,rows[0]);
          });

        },      
        function(err, results){        

          if (err) reject(err);

          gene_data['transcripts'] = results;

          ctx.set('Content-Type', 'application/json');
          ctx.set('Charset', 'utf-8')
          ctx.body = ctx.query.callback + '(' + JSON.stringify([gene_data]) +');';
          resolve();
        }
      );
    });
  });
});


router.get('/api/genes/', async (ctx) => {  

  var genesString = ctx.query.genes;
  var genes = genesString.split(",")

  var source = ctx.query.source;
  var species = ctx.query.species;
  var build = ctx.query.build;
  if (source == null || source == '') {
    source = 'gencode';
  } 
  var geneSqlString = "SELECT distinct * from genes where gene_name ";
  geneSqlString    += getGenesInClause(genes);
  geneSqlString    += " AND source = \""+source+"\"";
  if (species != null && species != "") {
    geneSqlString  += " AND species = \""+species+"\"";
  }
  if (build != null && build != "") {
    geneSqlString  += " AND build = \""+build+"\"";
  }

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(geneSqlString,function(err,rows){ 
      var gene_data = {};

      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          let gene_record = rows[i]
          gene_record.transcripts = []
          gene_data[gene_record.gene_name] = gene_record;    
        }
      } 
          
      var source = ctx.query.source; 
      if (source == null || source == '') {
        source = 'gencode';
      }
      var sqlString = "";
      sqlString =  "SELECT t.* ";
      sqlString += "FROM transcripts t "
      sqlString += "WHERE t.gene_name ";
      sqlString += getGenesInClause(genes);
      sqlString += " AND t.source = \""+source+"\"";
      if (species != null && species != "") {
        sqlString  += " AND t.species = \""+species+"\"";
      }
      if (build != null && build != "") {
        sqlString  += " AND t.build = \""+build+"\"";
      }        
      db.all(sqlString,function(err,transcriptRows){    

        if (err) reject(err);


        if (transcriptRows != null && transcriptRows.length > 0) {
          for (var i = 0; i < transcriptRows.length; i++) {
            var transcript = transcriptRows[i]; 
            transcript['features'] = JSON.parse(transcript['features']);
            gene_record = gene_data[transcript.gene_name]
            if (gene_record) {
              gene_record["transcripts"].push(transcript);
            } else {
              console.log("cannot find gene for transcript " + transcript.transcript_id + " " + transcript.gene_name)
            }  
          }
        }

        ctx.set('Content-Type', 'application/json');
        ctx.set('Charset', 'utf-8')
	ctx.body = JSON.stringify([gene_data]);
        resolve();
      });
    });
  });
});

router.get('/api/region/:region', async (ctx) => {  
  var chr = ctx.params.region.split(':')[0];
  var start = ctx.params.region.split(':')[1].split('-')[0];
  var end = ctx.params.region.split(':')[1].split('-')[1];
  var source = ctx.query.source; 
  var species = ctx.query.species;
  var build = ctx.query.build;
  
  // bound
  // 'outer'   (default) means start and end specified represent the outer-bounds.  i
  //           find all genes in the specified start and end region
  // 'inner'   means start and end specified represent a coordinate inside.
  //           in other words, find the gene that contains this start and end coordinate  
  var bound = ctx.query.bound;
  if (bound == null || bound == '') {
    bound = 'outer';
  }

  if (source == null || source == '') {
    source = 'gencode';
  } 
  var sqlString = "SELECT distinct * from genes where chr = '" + chr + "";
  if (bound == 'outer') {
    sqlString += "' and  (start between  " + start + " and " + end 
               + "        or end between " + start + " and " + end + ")";
  } else {
    sqlString += "' and  (start   <= " + start   
               + "        and end >= " + end + ")";
  }
  if (species != null && species != "") {
    sqlString  += " AND species = \""+species+"\"";
  }
  if (build != null && build != "") {
    sqlString  += " AND build = \""+build+"\"";
  }  
  if (source != null && source != "") {
    sqlString +=    " AND source = \""+source+"\"";         
  }
  
  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(sqlString, function(err, genes) {
      async.map(genes, 
        function(gene_data, outterDone) {                   
          var transcript_ids = JSON.parse(gene_data['transcripts']);
      
          async.map(transcript_ids,      
            function(id, done){      
              var sqlString = "SELECT * from transcripts t ";
              sqlString +=    "WHERE t.transcript_id=\""+id+"\" "
              if (source != null && source != "") {
		sqlString +=    " AND t.source = \""+source+"\""; 
	      }
              if (species != null && species != "") {
                sqlString  += " AND t.species = \""+species+"\"";
              }
              if (build != null && build != "") {
                sqlString  += " AND t.build = \""+build+"\"";
              }  
              db.all(sqlString,function(err,rows){          

                if (err) {
		  console.log("error: " + err);
		  reject(err);
		} 
                rows[0]['features'] = JSON.parse(rows[0]['features']);
                done(null,rows[0]);
              });
            },      
            function(err, results){        

              if (err) reject(err);

              gene_data['transcripts'] = results;            
              outterDone(null, gene_data);
            }
          );
        },
        function(err, results) {                

          if (err) reject(err);

          ctx.set('Content-Type', 'application/json');
          ctx.set('Charset', 'utf-8')
          ctx.body = JSON.stringify(results);
          resolve();
        }
      );
    }); 
  });
});



// v2 (cacheable) endpoints
router.get('/:gene', async (ctx) => {  
  var source = ctx.query.source;
  var species = ctx.query.species;
  var build = ctx.query.build;
  if (source == null || source == '') {
    source = 'gencode';
  } 
  var geneSqlString = "SELECT * from genes where gene_name like \""+ctx.params.gene+"\" ";
  geneSqlString    += " AND source = \""+source+"\"";
  if (species != null && species != "") {
    geneSqlString  += " AND species = \""+species+"\"";
  }
  if (build != null && build != "") {
    geneSqlString  += " AND build = \""+build+"\"";
  }

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(geneSqlString,function(err,rows){ 
      var gene_data = {};
      var transcript_ids = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          gene_data = rows[i];    
          if (gene_data.hasOwnProperty("transcripts") && gene_data.transcripts != null && gene_data.transcripts != "") {
            transcript_ids = transcript_ids.concat(JSON.parse(gene_data['transcripts']));
          }       
        }
      } 
          
      async.map(transcript_ids,      
        function(id, done){      
          var source = ctx.query.source; 
          if (source == null || source == '') {
            source = 'gencode';
          } 
          var sqlString = "";
          sqlString =   "SELECT t.* from transcripts t ";
          sqlString +=  "WHERE t.transcript_id=\""+id+"\" "
          sqlString +=  " AND t.source = \""+source+"\"";
          if (species != null && species != "") {
            sqlString  += " AND t.species = \""+species+"\"";
          }
          if (build != null && build != "") {
            sqlString  += " AND t.build = \""+build+"\"";
          }        
          db.all(sqlString,function(err,rows){    

            if (err) reject(err);

            if (rows != null && rows.length > 0) {
              rows[0]['features'] = JSON.parse(rows[0]['features']);
            } else {
              rows[0]['features'] = [];
            }   
            done(null,rows[0]);
          });

        },      
        function(err, results){        

          if (err) reject(err);

          gene_data['transcripts'] = results;

          ctx.set('Content-Type', 'application/json');
          ctx.set('Charset', 'utf-8')
          ctx.set('Cache-Control', 'public,max-age=84600')
          ctx.body = JSON.stringify([gene_data]);
          resolve();
        }
      );
    });
  });
});

/*
 * Return the transcript count by source and build
 * and return aliases for gene name.
 */
router.get('/lookupEntries/:genes', async (ctx) => {
  let geneWhereClause= ""
  let idx = 0;
  geneWhereClause = " g.gene_name in ("
  ctx.params.genes.split(",").forEach(function(geneName) {
    if (idx > 0) {
      geneWhereClause += ","
    }
    geneWhereClause += "'" + geneName + "'"
    idx++;
  })
  geneWhereClause +=  ")"

  let stmt = `
        SELECT
            g.gene_name,
            gs.gene_symbol,
            g.build,
            g.source,
            json_array_length(g.transcripts) as transcript_count,
            GROUP_CONCAT(ga.alias_symbol) AS aliases
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol and ga.alias_symbol != g.gene_name
        `

  stmt += " WHERE " + geneWhereClause
  stmt += " GROUP BY g.gene_name, gs.gene_symbol, g.build, g.source";

  console.log(stmt)
  const db = getDb();
  return new Promise((resolve, reject) => {
    db.all(stmt,function(err,rows){
      if (err) {
        console.log(err)
        reject(err);
      } else {
        console.log(rows)
        let gene_map = {}
	let gene_names = [];
        rows.forEach(function(row) {
          let gene_name        = row.gene_name
          let gene_symbol      = row.gene_symbol
          let build            = row.build
          let source           = row.source
          let transcript_count = row.transcript_count
          let aliases          = row.aliases

          if (gene_name == null || gene_name == "") {
            console.log("Warning, invalid gene name")
            console.log(row)
          }
          let gene = null;
          if (gene_map.hasOwnProperty(gene_name)) {
            gene = gene_map[gene_name]
          } else {
            gene ={   'gene_name': gene_name,
                      'GRCh37': {'gencode': 0, 'refseq': 0},
                      'GRCh38': {'gencode': 0, 'refseq': 0},
                   }
            gene_map[gene_name] = gene

            // Add the gene_symbol to the aliases (if the gene_symbol is
            // different than the gene name)
            if (gene_symbol && gene_symbol != gene_name) {
              if (aliases == null || aliases == "") {
                aliases = gene_symbol
              } else if (!aliases.hasOwnProperty(gene_symbol)) {
                aliases = gene_symbol + "," + aliases
              }
            }

            if (aliases && aliases.length > 0) {
              gene['aliases'] = aliases
            }

            if (gene_name && !gene_names.hasOwnProperty(gene_name)) {
              gene_names.push(gene_name)
            }
          }

          // Update the gene with the transcript count for the
          // row's build and source
          if (build && source && transcript_count) {
            gene[build][source] = transcript_count
          }
        }) // end of for loop over rows

        let genes = [];
        gene_names.forEach(function(geneName) {
         genes.push(gene_map[geneName])
        })
        var gene_data = {'genes': genes}
        ctx.set('Content-Type', 'application/json');
        ctx.set('Charset', 'utf-8')
        ctx.set('Cache-Control', 'public,max-age=84600')
        ctx.body = JSON.stringify(gene_data);
        resolve();
      } // end of else (not err)
    }) // end of db.all
  }) // end of new Promise
})

// Asynchronous lookup to support typeahead search based
// on all or part of gene name
router.get('/lookup/:gene', async (ctx) => {
  // searchAlias
  //   never  - Only search on gene names (known to gencode and refseq)
  //   always - Search on both gene names and aliases
  //   last   - Only search aliases if the searching on gene name returned no results
  // exactMatch
  //   true            - The gene name must match exactly (use = in WHERE clause)
  //   false (default) - The gene name starts with or equals the name provided
  var searchAlias = ctx.query.searchAlias;
  var exactMatch  = ctx.query.exactMatch && ctx.query.exactMatch == 'true' ? true : false;

  var stmt = "";

  if (searchAlias == 'always') {
    stmt = `SELECT distinct g.gene_name, ga.alias_symbol AS 'gene_alias'
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol and ga.alias_symbol != g.gene_name`
    if (exactMatch) {
      stmt += " WHERE g.gene_name     = " + "\"" + ctx.params.gene + "\""
      stmt += " OR    ga.alias_symbol = " + "\"" + ctx.params.gene + "\""
    } else {
      stmt += " WHERE g.gene_name     like " + "\"" + ctx.params.gene + "%" + "\""
      stmt += " OR    ga.alias_symbol like " + "\"" + ctx.params.gene + "%" + "\""
    }
    stmt += " GROUP BY g.gene_name"
  } else {
    stmt = "SELECT distinct g.gene_name from genes g "
    if (exactMatch) {
      stmt += " WHERE g.gene_name =    " + "\"" + ctx.params.gene + "\"";
    } else {
      stmt += " WHERE g.gene_name like " + "\"" + ctx.params.gene+"%" + "\"";
    }
  }

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(stmt,function(err,rows){
      if (err) {
        reject(err);
      } else {
        if (rows.length == 0 && searchAlias == 'last') {
          stmt = `SELECT distinct g.gene_name 'gene_name',
                                  ga.alias_symbol as 'gene_alias'
                  FROM      genes as g
                  LEFT JOIN gene_alias as ga on ga.gene_symbol = g.gene_symbol `
          if (exactMatch) {
            stmt += " WHERE ga.alias_symbol =    " + "\"" + ctx.params.gene + "\""
          } else {
            stmt += " WHERE ga.alias_symbol like " + "\"" + ctx.params.gene + "%" + "\""
          }
          db.all(stmt,function(err,rows){
            if (err) {
              reject(err);
            } else {
              var gene_data = {'genes': rows}
              ctx.set('Content-Type', 'application/json');
              ctx.set('Charset', 'utf-8')
              ctx.set('Cache-Control', 'public,max-age=84600')
              ctx.body = JSON.stringify(gene_data);
              resolve();
            }
          })
        } else {
          var gene_data = {'genes': rows};
          ctx.set('Content-Type', 'application/json');
          ctx.set('Charset', 'utf-8')
          ctx.set('Cache-Control', 'public,max-age=84600')
          ctx.body = JSON.stringify(gene_data);
          resolve();
        }
      }
    });
  });
});


module.exports = router;

