// Run with node tests/check_n8n_nodes.js. No SSH or n8n instance needed.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const workflow = JSON.parse(fs.readFileSync(path.join(__dirname, '../integrations/n8n/rnaseq.workflow.json')));
const code = name => workflow.nodes.find(n => n.name === name).parameters.jsCode;
const settings = new Function(code('Settings'))()[0].json;
assert.match(settings.startCommand, /--project/);
const lookup = () => ({first: () => ({json: settings})});
const input = json => ({first: () => ({json})});
const verify = new Function('$input', '$', code('Verify source version'));
assert.equal(verify(input({code:0,stdout:JSON.stringify({source_sha256:settings.expectedSource})}),lookup)[0].json.runId, settings.runId);
assert.throws(() => verify(input({code:0,stdout:'{"source_sha256":"old"}'}),lookup), /Wrong or modified/);
const parse = new Function('$input','$',code('Read status'));
for (const status of ['running','succeeded','failed','interrupted']) {
  const state = {project:settings.project,run_id:settings.runId,status};
  assert.equal(parse(input({code:0,stdout:JSON.stringify(state)}),lookup)[0].json.status,status);
}
for (const response of [{code:1,stderr:'error'},{code:0,stdout:'not json'},{code:0,stdout:'{}'}])
  assert.throws(() => parse(input(response),lookup));
const review = new Function('$input',code('Review results'));
assert.equal(review({all:()=>[{json:{mode:'run',reports:{recommendation:{route:'small_rna_plant'}}}}]})[0].json.route,'small_rna_plant');
assert.equal(review({all:()=>[{json:{mode:'raw_qc'}}]})[0].json.route,'independent_raw_qc');
console.log('PASS: n8n source gate, identity/error handling and route presentation');
const openSsh = JSON.parse(fs.readFileSync(path.join(__dirname, '../integrations/n8n/rnaseq.openssh.workflow.json')));
const sshCode = name => openSsh.nodes.find(n => n.name === name).parameters.jsCode;
const sshSettings = new Function(sshCode('Settings'))()[0].json;
assert.match(sshSettings.startCommand, /^\/usr\/bin\/ssh -o BatchMode=yes /);
assert.match(sshSettings.inspectCommand, /analysis-host/);
const sshLookup = () => ({first:()=>({json:sshSettings})});
const sshVerify = new Function('$input','$',sshCode('Verify source version'));
sshVerify(input({exitCode:0,stdout:JSON.stringify({source_sha256:sshSettings.expectedSource})}),sshLookup);
const sshParse = new Function('$input','$',sshCode('Read status'));
assert.equal(sshParse(input({exitCode:0,stdout:JSON.stringify({project:sshSettings.project,run_id:sshSettings.runId,status:'succeeded'})}),sshLookup)[0].json.status,'succeeded');
assert.throws(()=>sshParse(input({exitCode:255,stderr:'SSH failed'}),sshLookup));
console.log('PASS: OpenSSH command construction and exit-code handling');
