var sqlite3 = require('sqlite3').verbose();
const { dataPath } = require('./utils.js');
var db = new sqlite3.Database(dataPath('geneinfo/gene.iobio.db'));
var async = require('async');
const Router = require('koa-router');
const router = new Router();

router.get('/api/gene/:gene', async (ctx) => {  
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
          if (source == 'gencode') {
            sqlString =  "SELECT t.*, x.refseq_id as 'xref' from transcripts t ";
            sqlString += "LEFT OUTER JOIN xref_transcript x on x.gencode_id = t.transcript_id ";
          } else if (source == 'refseq') {
            sqlString =  "SELECT t.*, x.gencode_id as 'xref' from transcripts t ";
            sqlString += "LEFT OUTER JOIN xref_transcript x on x.refseq_id = t.transcript_id ";
          }
          sqlString +=    "WHERE t.transcript_id=\""+id+"\" "
          sqlString +=    "AND t.source = \""+source+"\"";
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
          //res.json([gene_data]);

          ctx.set('Content-Type', 'application/json');
          ctx.set('Charset', 'utf-8')
          ctx.body = ctx.query.callback + '(' + JSON.stringify([gene_data]) +');';
          resolve();
        }
      );
    });
  });
});

router.get('/api/region/:region', async (ctx) => {  
  var chr = ctx.params.region.split(':')[0].toLowerCase();
  var start = ctx.params.region.split(':')[1].split('-')[0];
  var end = ctx.params.region.split(':')[1].split('-')[1];
  var source = ctx.query.source; 
  var species = ctx.query.species;
  var build = ctx.query.build;
  if (source == null || source == '') {
    source = 'gencode';
  } 
  var sqlString = "SELECT * from genes where chr = '" + chr 
    + "' and  (start between " + start + " and " + end 
    + "        or end between " + start + " and " + end + ")";
  if (species != null && species != "") {
    sqlString  += " AND species = \""+species+"\"";
  }
  if (build != null && build != "") {
    sqlString  += " AND build = \""+build+"\"";
  }  

  return new Promise((resolve, reject) => {
    db.all(sqlString, function(err, genes) {  
      
      async.map(genes, 
        function(gene_data, outterDone) {                   
          var transcript_ids = JSON.parse(gene_data['transcripts']);
      
          async.map(transcript_ids,      
            function(id, done){      
              var sqlString = "SELECT * from transcripts t ";
              if (source == 'gencode') {
                sqlString +=    "LEFT OUTER JOIN xref_transcript x on x.gencode_id = t.transcript_id ";
              } else if (source == 'refseq') {
                sqlString +=    "LEFT OUTER JOIN xref_transcript x on x.refseq_id = t.transcript_id ";
              }
              sqlString +=    "WHERE t.transcript_id=\""+id+"\" "
              sqlString +=    "AND t.source = \""+source+"\"";          
              if (species != null && species != "") {
                sqlString  += " AND t.species = \""+species+"\"";
              }
              if (build != null && build != "") {
                sqlString  += " AND t.build = \""+build+"\"";
              }  
              db.all(sqlString,function(err,rows){          

                if (err) reject(err);

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

module.exports = router;
