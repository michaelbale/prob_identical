Shiny.addCustomMessageHandler("error",
  function(message) {
    alert(JSON.stringify(message));
  }
);