const fs = require('fs');
const html = fs.readFileSync('output_html_dump.html', 'utf8');
const scriptMatch = html.match(/<script>(.*?)<\/script>/s);
if(scriptMatch) {
    try {
        require('vm').Script(scriptMatch[1]);
        console.log("Valid JS");
    } catch(e) {
        console.log("INVALID JS:", e);
    }
}
