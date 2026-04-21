const puppeteer = require('puppeteer');
(async () => {
  const browser = await puppeteer.launch({args: ['--no-sandbox']});
  const page = await browser.newPage();
  
  page.on('console', msg => console.log('BROWSER LOG:', msg.text()));
  page.on('pageerror', err => console.log('BROWSER ERROR:', err.toString()));
  
  await page.goto('http://localhost:8501', {waitUntil: 'networkidle0'});
  
  // Wait a bit
  await new Promise(resolve => setTimeout(resolve, 5000));
  
  console.log("Dumping text:");
  const text = await page.evaluate(() => document.body.innerText);
  console.log(text.substring(0, 500));
  
  await browser.close();
})();
