const Router = require('koa-router');
var async = require('async');

let _db;
function getDb() {
  if (!_db) {
    const sqlite3 = require('sqlite3').verbose();
    const { dataPath } = require('./utils.js');
    _db = new sqlite3.Database(dataPath('hpo/hpo.db'));
  }
  return _db;
}

const router = new Router();

// TODO: commenting this out for now, since it's the same endpoint as the one
// below, so I'm pretty sure it's not doing anything. Need to verify which one
// we want then get rid of the other.
//
//router.get('/hpo/lookup', async (ctx) => {  
// 
//  var searchterm = ctx.query.term;
// 
//  var sqlString = "SELECT distinct HPO_DISEASE_NAME, SOURCE from hpo where HPO_DISEASE_NAME like \""+searchterm+"%\" order by HPO_DISEASE_NAME";
//  db.all(sqlString,function(err,rows){ 
//    var hpo_data= {};
//    var hpo_list=[];
//    if (rows != null && rows.length > 0) {
//      for (var i = 0; i < rows.length; i++) {
//        hpo_data= rows[i];           
//        hpo_list.push({id: hpo_data.SOURCE, label: hpo_data.HPO_DISEASE_NAME, value: hpo_data.HPO_DISEASE_NAME});
//      }
//    } 
//    ctx.set('Content-Type', 'application/json');
//    ctx.set('Charset', 'utf-8')
//    ctx.set("Access-Control-Allow-Origin","*")
//    ctx.set("Access-Control-Allow-Methods", "GET");
//    ctx.set("Access-Control-Allow-Headers", "Content-Type");
//    ctx.body = JSON.stringify(hpo_list);
//  });
//});

router.get('/hot/lookup', async (ctx) => {

  return new Promise((resolve, reject) => {
    var searchterm = ctx.query.term;

    var sqlString = "SELECT distinct disease_term from hot_disease_term where disease_term like \""+searchterm+"%\"";

    const db = getDb();
    
    db.all(sqlString,function(err,rows){
      var hot_data= {};
      var hot_list=[];
      if (rows != null && rows.length > 0) {
        for (var i = 0; i < rows.length; i++) {
          hot_data= rows[i];
          hot_list.push({id: hot_data.disease_term, label: hot_data.disease_term, value: hot_data.disease_term});
        }
      }
      ctx.set('Content-Type', 'application/json');
      ctx.set('Charset', 'utf-8')
      ctx.set("Access-Control-Allow-Origin","*")
      ctx.set("Access-Control-Allow-Methods", "GET");
      ctx.set("Access-Control-Allow-Headers", "Content-Type");
      ctx.body = JSON.stringify(hot_list);

      resolve();
    });
  });
});


module.exports = router;
