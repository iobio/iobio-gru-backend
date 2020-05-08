const Router = require('koa-router');
var sqlite3 = require('sqlite3').verbose();
const { dataPath } = require('./utils.js');
var db = new sqlite3.Database(dataPath('gene2pheno/hpo_gene_to_phenotype.db'));

const router = new Router();

router.get('/api/gene/:gene', async (ctx) => {
  var sqlString = "SELECT * from gene_to_phenotype where entrez_gene_symbol=\""+ctx.params.gene+"\" ";

  return new Promise((resolve, reject) => {
    db.all(sqlString,function(err,rows){ 

      if (err) reject(err);

      var phenotype_data = {};
      var phenotypes = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          phenotype_data = rows[i];           
          phenotypes.push(phenotype_data);
        }
      } 
      
      ctx.set('Content-Type', 'application/json');
      ctx.set('Charset', 'utf-8');
      ctx.body = ctx.query.callback + '(' + JSON.stringify(phenotypes) +');';

      resolve();
    });
  });
});


// v2 (cacheable) endpoints
router.get('/:gene', async (ctx) => {
  var sqlString = "SELECT * from gene_to_phenotype where entrez_gene_symbol=\""+ctx.params.gene+"\" ";

  return new Promise((resolve, reject) => {
    db.all(sqlString,function(err,rows){ 

      if (err) reject(err);

      var phenotype_data = {};
      var phenotypes = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          phenotype_data = rows[i];           
          phenotypes.push(phenotype_data);
        }
      } 
      
      ctx.set('Content-Type', 'application/json');
      ctx.set('Charset', 'utf-8');
      ctx.set('Cache-Control', 'public,max-age=84600')
      ctx.body = JSON.stringify(phenotypes);

      resolve();
    });
  });
});


module.exports = router;
