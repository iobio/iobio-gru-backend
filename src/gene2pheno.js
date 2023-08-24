const Router = require('koa-router');

let _db;
function getDb() {
  if (!_db) {
    const sqlite3 = require('sqlite3').verbose();
    const { dataPath } = require('./utils.js');
    _db = new sqlite3.Database(dataPath('gene2pheno/gene_to_phenotype.db'));
  }
  return _db;
}

function promiseGetGeneDisorders(db, gene_name) {
  var sqlString = "SELECT * from gene_to_disorder where gene_symbol=\""+gene_name+"\" ";
  return new Promise((resolve, reject) => {
    db.all(sqlString,function(err,rows){ 

      if (err) reject(err);

      var disorder_data = {};
      var disorders = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          disorder_data = rows[i];           
          disorders.push(disorder_data);
        }
      } 
      resolve(disorders)
    })
  })
}

const router = new Router();

router.get('/api/gene/:gene', async (ctx) => {
  var sqlString = "SELECT * from gene_to_phenotype where gene_symbol=\""+ctx.params.gene+"\" ";

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(sqlString,function(err,rows){ 

      if (err) reject(err);

      var phenotype_data = {};
      var phenotypes = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          phenotype_data = rows[i];  
          phenotype_data['entrez_gene_symbol'] = phenotype_data['gene_symbol']         
          phenotype_data['entrez_gene_id'] = phenotype_data['ncbi_gene_id']         
          phenotype_data['hpo_term_id'] = phenotype_data['hpo_id']         
          phenotype_data['hpo_term_name'] = phenotype_data['hpo_name']         
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
  var sqlString = "SELECT * from gene_to_phenotype where gene_symbol=\""+ctx.params.gene+"\" ";

  const db = getDb();

  return new Promise((resolve, reject) => {
    db.all(sqlString,function(err,rows){ 

      if (err) reject(err);

      var phenotype_data = {};
      var phenotypes = [];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          phenotype_data = rows[i];    
          phenotype_data['entrez_gene_symbol'] = phenotype_data['gene_symbol']         
          phenotype_data['entrez_gene_id'] = phenotype_data['ncbi_gene_id']         
          phenotype_data['hpo_term_id'] = phenotype_data['hpo_id']         
          phenotype_data['hpo_term_name'] = phenotype_data['hpo_name']         
                 
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


// v3 In addition to returning associated phenotypes, this endpoint also returns 
// associated disorders for a gene.
router.get('/associations/:gene', async (ctx) => {
  var sqlString = "SELECT * from gene_to_phenotype where gene_symbol=\""+ctx.params.gene+"\" ";

  const db = getDb();

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

      promiseGetGeneDisorders(db, ctx.params.gene)
      .then(function(disorders) {

        let associations = {'gene': ctx.params.gene, 
                            'phenotypes': phenotypes, 
                            'disorders': disorders}
      
        ctx.set('Content-Type', 'application/json');
        ctx.set('Charset', 'utf-8');
        ctx.set('Cache-Control', 'public,max-age=84600')
        ctx.body = JSON.stringify(associations);

        resolve();

      })
      .catch(function(error) {
        reject(error)
      })


    });
  });
});



module.exports = router;
