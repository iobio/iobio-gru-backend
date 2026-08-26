import sqlite3 from 'sqlite3';
import { dataPath } from './utils.js';

let _db;
function getDb() {
  if (!_db) {
    const sqlite3Verbose = sqlite3.verbose();
    _db = new sqlite3Verbose.Database(dataPath('geneinfo/gene.iobio.db'));
  }
  return _db;
}

function dbAll(sql) {
  const db = getDb();
  return new Promise((resolve, reject) => {
    db.all(sql, (err, rows) => {
      if (err) reject(err);
      else resolve(rows);
    });
  });
}

function jsonResponse(body, opts = {}) {
  const headers = {
    'Content-Type': 'application/json',
    'Charset': 'utf-8',
  };

  if (opts.cache) {
    headers['Cache-Control'] = 'public,max-age=84600';
  }

  const init = { headers };
  if (opts.status) {
    init.status = opts.status;
  }

  return new Response(body, init);
}

function badRequest(message) {
  return jsonResponse(JSON.stringify({ error: message }), { status: 400 });
}

function getGenesInClause(genes) {
  let sqlString = ' in (';
  let firstTime = true;

  genes.forEach(function(gene) {
    if (firstTime) {
      firstTime = false;
    } else {
      sqlString += ',';
    }
    sqlString += '"' + gene + '" ';
  });

  sqlString += ')';
  return sqlString;
}

function getPathPart(pathParts, idx) {
  try {
    return decodeURIComponent(pathParts[idx] || '');
  } catch (err) {
    return null;
  }
}

function createGeneInfoHandler(opt) {
  return async (req, url = new URL(req.url)) => {
    try {
      const params = Object.fromEntries(url.searchParams);

      const path = url.pathname.slice(opt.prefix.length);
      const pathParts = path.split('/');
      const command = pathParts[1];

      let response;
      if (command === 'api') {
        const apiCommand = pathParts[2];

        if (apiCommand === 'gene') {
          response = apiGene(getPathPart(pathParts, 3), params);
        }
        else if (apiCommand === 'genes') {
          response = apiGenes(params);
        }
        else if (apiCommand === 'region') {
          response = apiRegion(getPathPart(pathParts, 3), params);
        }
        else if (apiCommand === 'lookupGenes') {
          response = apiLookupGenes(params);
        }
        else {
          response = new Response('Not found', { status: 404 });
        }
      }
      else if (command === 'lookupEntries') {
        response = lookupEntries(getPathPart(pathParts, 2));
      }
      else if (command === 'lookup') {
        response = lookup(getPathPart(pathParts, 2), params);
      }
      else {
        response = geneinfo(getPathPart(pathParts, 1), params);
      }

      return await response;
    } catch (err) {
      console.error('Unhandled geneinfo error:', err);
      return jsonResponse(JSON.stringify({ error: 'Internal server error' }), { status: 500 });
    }
  };
}

// Legacy /api/* endpoints retained for clients that still call the pre-v2 API.
async function apiGene(gene, params) {
  if (gene == null || gene == '') {
    return badRequest('Missing required path parameter: gene');
  }

  let source = params.source;
  const { species, build } = params;

  if (source == null || source == '') {
    source = 'gencode';
  }

  let geneSqlString = 'SELECT g.*, ';
  geneSqlString += 'gs.gene_symbol, ';
  geneSqlString += 'GROUP_CONCAT(distinct ga.alias_symbol) AS aliases ';
  geneSqlString += 'FROM genes g ';
  geneSqlString += 'LEFT OUTER JOIN gene_symbol gs on gs.gene_symbol = g.gene_symbol ';
  geneSqlString += 'LEFT OUTER JOIN gene_alias  ga on ga.gene_symbol = g.gene_symbol and ga.alias_symbol != g.gene_name ';
  geneSqlString += 'WHERE gene_name like "' + gene + '" ';
  geneSqlString += ' AND source = "' + source + '"';
  if (species != null && species != '') {
    geneSqlString  += ' AND species = "' + species + '"';
  }
  if (build != null && build != '') {
    geneSqlString  += ' AND build = "' + build + '"';
  }

  const rows = await dbAll(geneSqlString);
  let gene_data = {};
  let transcript_ids = [];

  if (rows != null && rows.length > 0) {
    for (let i = 0; i < rows.length; i++) {
      gene_data = rows[i];
      if (gene_data.hasOwnProperty('transcripts') && gene_data.transcripts != null && gene_data.transcripts != '') {
        transcript_ids = transcript_ids.concat(JSON.parse(gene_data.transcripts));
      }
    }
  }

  const transcripts = await Promise.all(transcript_ids.map(async (id) => {
    let sqlString = '';
    sqlString =  'SELECT t.* ';
    sqlString += 'from transcripts t ';
    sqlString += 'WHERE t.transcript_id="' + id + '" ';
    sqlString += ' AND t.source = "' + source + '"';
    if (species != null && species != '') {
      sqlString  += ' AND t.species = "' + species + '"';
    }
    if (build != null && build != '') {
      sqlString  += ' AND t.build = "' + build + '"';
    }

    const transcriptRows = await dbAll(sqlString);
    if (transcriptRows != null && transcriptRows.length > 0) {
      transcriptRows[0].features = JSON.parse(transcriptRows[0].features);
      return transcriptRows[0];
    }

    return [];
  }));

  gene_data.transcripts = transcripts;

  return jsonResponse(params.callback + '(' + JSON.stringify([gene_data]) + ');');
}

async function apiGenes(params) {
  const genesString = params.genes;
  if (genesString == null || genesString == '') {
    return badRequest('Missing required query parameter: genes');
  }

  const genes = genesString.split(',');

  let source = params.source;
  const { species, build } = params;

  if (source == null || source == '') {
    source = 'gencode';
  }

  let geneSqlString = 'SELECT distinct * from genes where gene_name ';
  geneSqlString    += getGenesInClause(genes);
  geneSqlString    += ' AND source = "' + source + '"';
  if (species != null && species != '') {
    geneSqlString  += ' AND species = "' + species + '"';
  }
  if (build != null && build != '') {
    geneSqlString  += ' AND build = "' + build + '"';
  }

  const rows = await dbAll(geneSqlString);
  const gene_data = {};

  if (rows != null && rows.length > 0) {
    for (let i = 0; i < rows.length; i++) {
      const gene_record = rows[i];
      gene_record.transcripts = [];
      gene_data[gene_record.gene_name] = gene_record;
    }
  }

  let sqlString = '';
  sqlString =  'SELECT t.* ';
  sqlString += 'FROM transcripts t ';
  sqlString += 'WHERE t.gene_name ';
  sqlString += getGenesInClause(genes);
  sqlString += ' AND t.source = "' + source + '"';
  if (species != null && species != '') {
    sqlString  += ' AND t.species = "' + species + '"';
  }
  if (build != null && build != '') {
    sqlString  += ' AND t.build = "' + build + '"';
  }

  const transcriptRows = await dbAll(sqlString);

  if (transcriptRows != null && transcriptRows.length > 0) {
    for (let i = 0; i < transcriptRows.length; i++) {
      const transcript = transcriptRows[i];
      transcript.features = JSON.parse(transcript.features);
      const gene_record = gene_data[transcript.gene_name];
      if (gene_record) {
        gene_record.transcripts.push(transcript);
      } else {
        console.log('cannot find gene for transcript ' + transcript.transcript_id + ' ' + transcript.gene_name);
      }
    }
  }

  return jsonResponse(JSON.stringify([gene_data]));
}

async function apiRegion(region, params) {
  if (region == null || region == '') {
    return badRequest('Missing required path parameter: region');
  }

  const regionParts = region.split(':');
  const coordParts = regionParts[1] ? regionParts[1].split('-') : [];
  if (regionParts.length < 2 || coordParts.length < 2 || !regionParts[0] || !coordParts[0] || !coordParts[1]) {
    return badRequest('Invalid region. Expected chr:start-end');
  }

  const chr = regionParts[0];
  const start = coordParts[0];
  const end = coordParts[1];
  if (!/^\d+$/.test(start) || !/^\d+$/.test(end)) {
    return badRequest('Invalid region coordinates. Expected numeric start and end');
  }
  let source = params.source;
  const { species, build } = params;

  // bound
  // 'outer'   (default) means start and end specified represent the outer-bounds.  i
  //           find all genes in the specified start and end region
  // 'inner'   means start and end specified represent a coordinate inside.
  //           in other words, find the gene that contains this start and end coordinate
  let bound = params.bound;
  if (bound == null || bound == '') {
    bound = 'outer';
  }

  if (source == null || source == '') {
    source = 'gencode';
  }

  let sqlString = "SELECT distinct * from genes where chr = '" + chr + '';
  if (bound == 'outer') {
    sqlString += "' and  (start between  " + start + ' and ' + end
               + '        or end between ' + start + ' and ' + end + ')';
  } else {
    sqlString += "' and  (start   <= " + start
               + '        and end >= ' + end + ')';
  }
  if (species != null && species != '') {
    sqlString  += ' AND species = "' + species + '"';
  }
  if (build != null && build != '') {
    sqlString  += ' AND build = "' + build + '"';
  }
  if (source != null && source != '') {
    sqlString += ' AND source = "' + source + '"';
  }

  const genes = await dbAll(sqlString);
  const results = await Promise.all((genes || []).map(async (gene_data) => {
    const transcript_ids = JSON.parse(gene_data.transcripts);

    const transcripts = await Promise.all(transcript_ids.map(async (id) => {
      let transcriptSqlString = 'SELECT * from transcripts t ';
      transcriptSqlString += 'WHERE t.transcript_id="' + id + '" ';
      if (source != null && source != '') {
        transcriptSqlString += ' AND t.source = "' + source + '"';
      }
      if (species != null && species != '') {
        transcriptSqlString  += ' AND t.species = "' + species + '"';
      }
      if (build != null && build != '') {
        transcriptSqlString  += ' AND t.build = "' + build + '"';
      }

      const rows = await dbAll(transcriptSqlString);
      if (rows != null && rows.length > 0) {
        rows[0].features = JSON.parse(rows[0].features);
        return rows[0];
      }

      return [];
    }));

    gene_data.transcripts = transcripts;
    return gene_data;
  }));

  return jsonResponse(JSON.stringify(results));
}

// v2 (cacheable) endpoints
async function geneinfo(gene, params) {
  if (gene == null || gene == '') {
    return badRequest('Missing required path parameter: gene');
  }

  const source = params.source || 'gencode';
  const { species, build } = params;

  let geneSqlString = `SELECT * from genes where gene_name like "${gene}" `;
  geneSqlString += ` AND source = "${source}"`;
  if (species) geneSqlString += ` AND species = "${species}"`;
  if (build) geneSqlString += ` AND build = "${build}"`;

  const geneRows = await dbAll(geneSqlString);

  let gene_data = {};
  let transcript_ids = [];

  if (geneRows && geneRows.length > 0) {
    for (const row of geneRows) {
      gene_data = row;
      if (row.transcripts) {
        transcript_ids = transcript_ids.concat(
          JSON.parse(row.transcripts)
        );
      }
    }
  }

  // Fetch transcripts in parallel
  const transcripts = await Promise.all(
    transcript_ids.map(async (id) => {
      let sqlString = `
        SELECT t.* from transcripts t
        WHERE t.transcript_id="${id}"
        AND t.source="${source}"
      `;
      if (species) sqlString += ` AND t.species="${species}"`;
      if (build) sqlString += ` AND t.build="${build}"`;

      const rows = await dbAll(sqlString);

      if (rows && rows.length > 0) {
        rows[0].features = JSON.parse(rows[0].features);
        return rows[0];
      }

      return {};
    })
  );

  gene_data.transcripts = transcripts;

  return jsonResponse(JSON.stringify([gene_data]), { cache: true });
}

/*
 * Return the transcript count by source and build
 * and return aliases for gene name.
 */
function lookupEntries(genes) {
  if (genes == null || genes == '') {
    return badRequest('Missing required path parameter: genes');
  }

  let geneWhereClause= ""
  let idx = 0;
  geneWhereClause = " g.gene_name in ("
  genes.split(",").forEach(function(geneName) {
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
            g.chr,
            g.start,
            g.end,
            json(g.transcripts) as transcript_ids,
            GROUP_CONCAT(ga.alias_symbol) AS aliases
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol and ga.alias_symbol != g.gene_name
        `

  stmt += " WHERE " + geneWhereClause
  stmt += " GROUP BY g.gene_name, gs.gene_symbol, g.build, g.source, g.chr, g.start, g.end";
  const db = getDb();
  return new Promise((resolve, reject) => {
    db.all(stmt,function(err,rows){
      if (err) {
        console.log(err)
        reject(err);
      } else {
        let gene_map = {}
	let gene_names = [];
        rows.forEach(function(row) {
          let gene_name        = row.gene_name
          let gene_symbol      = row.gene_symbol
          let build            = row.build
          let source           = row.source
          let chr              = row.chr
          let start            = row.start
          let end              = row.end
          let transcript_ids   = row.transcript_ids
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
                      'gene_coord': {
                        'GRCh37': {'gencode': {}, 'refseq': {}},
                        'GRCh38': {'gencode': {}, 'refseq': {}},
                      }
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

          let valid_transcript_ids  = [];
          if (transcript_ids && transcript_ids != '[null]') {
            valid_transcript_ids = JSON.parse(transcript_ids)
          }


          // Update the gene with the transcript count for the
          // row's build and source
          if (build && source) {
            gene[build][source] = valid_transcript_ids.length
          }
          // Update the gene coords
          if (build && source) {
            let coord = {"chr": chr, "start": start, "end": end};
            gene.gene_coord[build][source] = coord;
          }
        }) // end of for loop over rows

        let genes = [];
        gene_names.forEach(function(geneName) {
         genes.push(gene_map[geneName])
        })
        var gene_data = {'genes': genes}

        resolve(jsonResponse(JSON.stringify(gene_data), { cache: true }));
      } // end of else (not err)
    }) // end of db.all
  }) // end of new Promise
}

// Asynchronous lookup to support typeahead search based
// on all or part of gene name
function lookup(gene, params) {
  if (gene == null || gene == '') {
    return badRequest('Missing required path parameter: gene');
  }

  // searchAlias
  //   never  - Only search on gene names (known to gencode and refseq)
  //   always - Search on both gene names and aliases
  //   last   - Only search aliases if the searching on gene name returned no results
  // exactMatch
  //   true            - The gene name must match exactly (use = in WHERE clause)
  //   false (default) - The gene name starts with or equals the name provided
  var searchAlias = params.searchAlias;
  var exactMatch  = params.exactMatch && params.exactMatch == 'true' ? true : false;

  var stmt = "";

  if (searchAlias == 'always') {
    stmt = `SELECT distinct g.gene_name, ga.alias_symbol AS 'gene_alias'
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol and ga.alias_symbol != g.gene_name`
    if (exactMatch) {
      stmt += " WHERE g.gene_name     = " + "\"" + gene + "\""
      stmt += " OR    ga.alias_symbol = " + "\"" + gene + "\""
    } else {
      stmt += " WHERE g.gene_name     like " + "\"" + gene + "%" + "\""
      stmt += " OR    ga.alias_symbol like " + "\"" + gene + "%" + "\""
    }
    stmt += " GROUP BY g.gene_name"
  } else {
    stmt = "SELECT distinct g.gene_name from genes g "
    if (exactMatch) {
      stmt += " WHERE g.gene_name =    " + "\"" + gene + "\"";
    } else {
      stmt += " WHERE g.gene_name like " + "\"" + gene+"%" + "\"";
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
            stmt += " WHERE ga.alias_symbol =    " + "\"" + gene + "\""
          } else {
            stmt += " WHERE ga.alias_symbol like " + "\"" + gene + "%" + "\""
          }
          db.all(stmt,function(err,rows){
            if (err) {
              reject(err);
            } else {
              var gene_data = {'genes': rows}
              resolve(jsonResponse(JSON.stringify(gene_data), { cache: true }));
            }
          })
        } else {
          var gene_data = {'genes': rows};
          resolve(jsonResponse(JSON.stringify(gene_data), { cache: true }));
        }
      }
    });
  });
}

// Lookup multiple genes and return an array of the genes
// found, either as the gene name or an alias. Use like
// in the where clause to perform a case-insensitive search.
async function apiLookupGenes(params) {
  // genes    - The comma separated list of genes
  // searchAlias
  //   never  - Only search on gene names (known to gencode and refseq)
  //   always - Search on both gene names and aliases
  //   last   - Only search aliases if the searching on gene name returned no results
  const genes       = params.genes;
  if (genes == null || genes == '') {
    return badRequest('Missing required query parameter: genes');
  }

  const searchAlias = params.searchAlias;
  let stmt = '';

  if (searchAlias == 'always') {
    stmt = `SELECT distinct g.gene_name, ga.alias_symbol AS 'gene_alias'
        FROM genes g
        LEFT JOIN gene_symbol gs
          ON g.gene_symbol = gs.gene_symbol
        LEFT JOIN gene_alias ga
          ON gs.gene_symbol = ga.gene_symbol`

    let firstTime = true;
    genes.split(',').forEach(function(geneName) {
      stmt += firstTime ? ' WHERE ' : ' OR '
      stmt += ' g.gene_name     like ' + '"' + geneName + '"'
      stmt += ' OR    ga.alias_symbol like ' + '"' + geneName + '"'
      firstTime = false;
    })
    stmt += ' GROUP BY g.gene_name'
  } else {
    stmt = 'SELECT distinct g.gene_name from genes g '
    let firstTime = true;
    genes.split(',').forEach(function(geneName) {
      stmt += firstTime ? ' WHERE ' : ' OR '
      stmt += ' g.gene_name     like ' + '"' + geneName + '"'
      firstTime = false;
    })
  }

  const rows = await dbAll(stmt);

  // Map gene name parameter to result row
  const inputGeneMap = {};
  genes.split(',').forEach(function(inputGeneName) {
    const inputGeneNameUC = inputGeneName.toUpperCase();
    inputGeneMap[inputGeneNameUC] = false;
  })

  rows.forEach(function(row) {
    if (!row.gene_name) {
      return;
    }

    let matchedRow = inputGeneMap[row.gene_name.toUpperCase()]
    if (matchedRow != null && matchedRow == false) {
      const result = {'match': true, 'gene_name': row.gene_name};
      inputGeneMap[row.gene_name.toUpperCase()] = result;
    } else if (row.gene_alias != null) {
      matchedRow = inputGeneMap[row.gene_alias.toUpperCase()]
      if (matchedRow != null && matchedRow == false) {
        const result = {'match': true, 'gene_alias': row.gene_alias};
        inputGeneMap[row.gene_alias.toUpperCase()] = result;
      }
    }
  })

  const genesNoMatch = [];
  genes.split(',').forEach(function(inputGene) {
    const matchedRow = inputGeneMap[inputGene.toUpperCase()];
    if (!matchedRow) {
      genesNoMatch.push(inputGene)
    }
  })

  if (genesNoMatch.length > 0 && searchAlias == 'last') {
    stmt = `SELECT distinct g.gene_name 'gene_name',
                            ga.alias_symbol as 'gene_alias'
            FROM      genes as g
            LEFT JOIN gene_alias as ga on ga.gene_symbol = g.gene_symbol `

    let firstTime = true;
    genesNoMatch.forEach(function(geneName) {
      stmt += firstTime ? ' WHERE ' : ' OR '
      stmt += ' ga.alias_symbol like ' + '"' + geneName  + '"'
      firstTime = false;
    })

    const aliasRows = await dbAll(stmt);

    // Match the input gene name to the result row for gene alias
    aliasRows.forEach(function(aliasRow) {
      if (aliasRow.gene_alias == null) {
        return;
      }

      const matchedRow = inputGeneMap[aliasRow.gene_alias.toUpperCase()]
      if (matchedRow != null && matchedRow == false) {
        const result = {'match': true, 'gene_alias': aliasRow.gene_alias}
        inputGeneMap[aliasRow.gene_alias.toUpperCase()] = result;
      }
    })

    // Now create a result set for every gene name in the parameter
    const results = []
    genes.split(',').forEach(function(inputGene) {
      let resultRow = inputGeneMap[inputGene.toUpperCase()];
      if (resultRow) {
        resultRow.input_gene_name = inputGene;
      } else {
        resultRow = {'input_gene_name': inputGene, 'gene_name': '', 'gene_alias': '', 'match': ''}
      }
      results.push(resultRow)
    })

    return jsonResponse(JSON.stringify({'genes': results}), { cache: true });
  } else {
    // Now create a result set for every gene name in the parameter
    const results = []
    genes.split(',').forEach(function(inputGene) {
      let resultRow = inputGeneMap[inputGene.toUpperCase()];
      if (resultRow) {
        resultRow.input_gene_name = inputGene;
      } else {
        resultRow = {'match': false, 'input_gene_name': inputGene};
      }
      results.push(resultRow)
    })

    return jsonResponse(JSON.stringify({'genes': results}), { cache: true });
  }
}

export {
  createGeneInfoHandler,
};
