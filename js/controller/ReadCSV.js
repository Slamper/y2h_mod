xinet.Controller.prototype.readCSV = function (csvContents) {
    var rows = d3.csv.parseRows(csvContents);
    var $j = jQuery.noConflict();
    var headers = rows[0];
    console.log(headers.toString());

    var uniProt = headers.indexOf('UniProt');
    var iId = headers.indexOf('ORF');
    var iCat = headers.indexOf('Category');

    //missing UniProt
    if (uniProt === -1) {
        alert("Failed to read column 'UniProt' from CSV file");
        return;
    }
    if (iId === -1) {
        alert("Failed to read column 'ID' from CSV file");
        return;
    }
    if (iCat === -1) {
        alert("Failed to read column 'Category' from CSV file");
        return;
    }

    // if a Fasta file has been added then this.proteins will not be empty
    var countRows = rows.length;
    if (this.proteins.keys().length === 0) {
        //No protein data. We will need to look up accession numbers to get sequences.

        //We are likely going to encounter things like proteins with
        //differnt ids/names but the same accession number.
        var accLookupMap = d3.map();
        //The following server is not returning results for protein isoforms.
        var server_url = 'http://www.ebi.ac.uk/das-srv/uniprot/das/uniprot/';
        var client = JSDAS.Simple.getClient(server_url);
        addProteins(this);
        initProteins(this);
    }

    function addProteins(xlv) {
        for (var row = 1; row < countRows; row++) {
            var protsUni = rows[row][uniProt];
            var protsName = rows[row][iId];
            var cat = rows[row][iCat];
            var id = protsName,
                acc = protsUni,
                name = protsName;
            if (!xlv.proteins.has(id)) {
                var protein = new Protein(id, xlv, acc, name, cat);
                xlv.proteins.set(id, protein);
                var accLookupEntry = accLookupMap.get(acc);
                if (typeof accLookupEntry === "undefined") {
                    accLookupMap.set(acc, [id]);
                } else {
                    accLookupEntry.push(id);
                }
            }
        }
    }

    function initProteins(xlv) {
        // This function will be executed in case of error
        var error_response = function () {
            alert('No FASTA file and DAS sequence look up failed.');
        };
        // This function inits the protein
        var response = function (res) {
            //this.message(res);
            var acc = res.SEQUENCE[0].id.trim();
            var seq = res.SEQUENCE[0].textContent.trim();
            var label = res.SEQUENCE[0].label.trim();
            var pids = accLookupMap.get(acc);
            for (var i = 0; i < pids.length; i++) {
                var prot = xlv.proteins.get(pids[i]);
                prot.initProtein(seq, label);
            }
            accLookupMap.remove(acc);
            xlv.message('Waiting on DAS response (sequence) for:<br/>' + accLookupMap.keys().toString());
            if (accLookupMap.keys().length === 0) {
                xlv.message('All sequences downloaded from DAS');
                addCSVLinks(xlv);
            }
        };
        var keys = accLookupMap.keys();
        var accCount = keys.length;
        for (var p = 0; p < accCount; p++) {
            var accession = keys[p];
            if (accession !== "") {
                //Asking the client to retrieve the sequence
                xiNET_Storage.getSequence(accession, function (ident, seq) {
                    //this.message(res);
                    var acc = ident;
                    var seq = seq.trim();
                    var label = "";
                    var pids = accLookupMap.get(acc);
                    for (var i = 0; i < pids.length; i++) {
                        var prot = xlv.proteins.get(pids[i]);
                        prot.initProtein(seq, label);
                    }
                    accLookupMap.remove(acc);
                    xlv.message('Waiting on DAS response (sequence) for:<br/>' + accLookupMap.keys().toString());
                    if (accLookupMap.keys().length === 0) {
                        xlv.message('All sequences downloaded from DAS');
                        addCSVLinks(xlv);
                    }
                    }
                );
                //client.sequence({
                //    segment: accession
                //}, response, error_response);
            }
            else {


                accLookupMap.remove(accession);
            }
        }
    }

    function addCSVLinks(xlv) {
        var prot1, prot2, id, score;
        for (var row = 1; row < countRows; row++) {
            for (var column = 3; column < headers.length; column++) {
                prot1 = rows[row][iId];
                prot2 = headers[column];
                score = rows[row][column];
                if (score != "0" && score != "") {
                    xlv.addMatch(prot1, "match-prot1",
                        prot2, "",
                        "match-" + row + "-" + column, score);
                }
            }
        }
        //~ }
        var protCount = xlv.proteins.values().length;
        var prots = xlv.proteins.values();
        /*for (var p = 0; p < protCount; p++) {
         var prot = prots[p];
         if (*/
        /*prot.name.indexOf("DECOY_") !== -1 &&*/
        /* prot.proteinLinks.keys().length === 0) {
         xlv.proteins.remove(prot.id);
         }
         }     */

        xlv.init();
        var keys = xlv.proteins.keys();
        $j("#search").autocomplete({
            source: xlv.proteins.keys(),
            appendTo: "#controls"
        });
        if (typeof initSlider === "function") {
            initSlider();
        }
        new xinet.DASUtil(xlv);
    }
};
//~ 
//~ xinet.Controller.prototype.readXQuest = function(csvContents) {
//~ var rows = d3.csv.parse(csvContents);
//~ //    var headers = rows[0];//first row is headers
//~ //    var iProt1 = headers.indexOf('protein1');
//~ //    var iRes1 = headers.indexOf('residue1');
//~ //    var iProt2 = headers.indexOf('protein2');
//~ //    var iRes2 = headers.indexOf('residue2');
//~ //    var iDescription = headers.indexOf('description');
//~ var countRows = rows.length;
//~ var prot1, prot2;
//~ for (var row = 0; row < countRows; row++) {
//~ prot1 = rows[row]['Protein1'].trim();
//~ prot2 = rows[row]['Protein2'].trim();
//~ if (prot1.toLowerCase().indexOf("reverse") === -1 && prot2.toLowerCase().indexOf("reverse") === -1
//~ && prot1.toLowerCase().indexOf("decoy") === -1 && prot2.toLowerCase().indexOf("decoy") === -1) {
//~ xlv.addMatch(prot1, rows[row]['AbsPos1'], prot2, rows[row]['AbsPos2'],
//~ rows[row]['Id'], rows[row]['ld-Score']);
//~ }
//~ }
//~ xlv.init();
//~ if (typeof initSlider === "function"){
//~ initSlider();
//~ }
//~ };
